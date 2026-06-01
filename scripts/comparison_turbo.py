"""Turbo equalization BER comparison: Ungerboeck M-BCJR vs Paper [14] vs Paper [26].

Uses (7,5) rate-1/2 convolutional code with iterative equalization/decoding.
All algorithms run under identical FTN channel conditions.

Curves:
  - No-ISI baseline (AWGN)
  - Ungerboeck M-BCJR (M=4, M=8)
  - Paper [14] Smoothed Backup M-BCJR (M=4, M=8)
  - Paper [26] Shortened BCJR (nu=2, nu=3)

Output: results/comparison/turbo_tau{05,035}_ber.csv + PNG
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ftn.baselines.channel_shortening import (
    apply_shortening_filter_fft,
    compute_shortened_params_from_v,
)
from ftn.baselines.forney_model import forney_channel, min_phase_from_pulse
from ftn.baselines.prlja_mbcjr import prlja_backup_mbcjr_bpsk
from ftn.baselines.shortened_bcjr import shortened_bcjr_bpsk
from ftn.channel import ftn_awgn_channel
from ftn.coding.conv import conv_bcjr_decode, conv_encode_75
from ftn.equalizers.ungerboeck_mbcjr import ungerboeck_mbcjr_bpsk
from ftn.modulation import bpsk_modulate
from ftn.pulse import compute_g, generate_rrc


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_interleaver(n, seed):
    return np.random.default_rng(seed).permutation(n)


def _noisi_frame(info_bits, es_n0_db, rng):
    n0 = 1.0 / (10.0 ** (es_n0_db / 10.0))
    encoded = conv_encode_75(info_bits)
    symbols = bpsk_modulate(encoded)
    noise = rng.standard_normal(symbols.size) * np.sqrt(n0 / 2.0)
    y = symbols + noise
    code_llr = 4.0 * y / n0
    result = conv_bcjr_decode(code_llr)
    hard_info = (result.info_llr < 0.0).astype(np.uint8)
    errors = int(np.sum(hard_info != info_bits))
    return errors, int(info_bits.size)


# ---------------------------------------------------------------------------
# Turbo frames for each equalizer
# ---------------------------------------------------------------------------

def _ungerboeck_turbo_frame(info_bits, g, n0, isi_len, m_states, future_len,
                            turbo_iters, llr_clip, interleaver, rng):
    encoded = conv_encode_75(info_bits)
    inv_perm = np.argsort(interleaver)
    int_encoded = encoded[interleaver]
    symbols = bpsk_modulate(int_encoded)
    y = ftn_awgn_channel(symbols, g, n0=n0, rng=rng)

    detector_prior = np.zeros(y.size, dtype=float)
    initial_state = tuple(1.0 for _ in range(isi_len))

    for _ in range(turbo_iters + 1):
        det = ungerboeck_mbcjr_bpsk(
            y, g, n0=n0, la=detector_prior,
            isi_len=isi_len, m_states=m_states,
            future_len=future_len, initial_state=initial_state,
        )
        det_llr = np.clip(det.llr, -llr_clip, llr_clip)
        det_ext = np.clip(det_llr - detector_prior, -llr_clip, llr_clip)
        det_ext_deint = det_ext[inv_perm]

        dec = conv_bcjr_decode(det_ext_deint)
        info_llr = np.clip(dec.info_llr, -llr_clip, llr_clip)
        code_llr = np.clip(dec.code_llr, -llr_clip, llr_clip)
        dec_ext = np.clip(code_llr - det_ext_deint, -llr_clip, llr_clip)
        detector_prior = dec_ext[interleaver]

    hard_info = (info_llr < 0.0).astype(np.uint8)
    errors = int(np.sum(hard_info != info_bits))
    return errors, int(info_bits.size)


def _prlja_turbo_frame(info_bits, v, n0, M, turbo_iters, llr_clip,
                       interleaver, rng):
    encoded = conv_encode_75(info_bits)
    inv_perm = np.argsort(interleaver)
    int_encoded = encoded[interleaver]
    symbols = bpsk_modulate(int_encoded)
    # Forney model: y = conv(symbols, v) + white_noise
    y = forney_channel(symbols, v, n0=n0, rng=rng)

    detector_prior = np.zeros(y.size, dtype=float)

    for _ in range(turbo_iters + 1):
        det = prlja_backup_mbcjr_bpsk(
            y, v, n0=n0, M=M, M_B=2, smooth=True, la=detector_prior,
        )
        det_llr = np.clip(det.llr, -llr_clip, llr_clip)
        det_ext = np.clip(det_llr - detector_prior, -llr_clip, llr_clip)
        det_ext_deint = det_ext[inv_perm]

        dec = conv_bcjr_decode(det_ext_deint)
        info_llr = np.clip(dec.info_llr, -llr_clip, llr_clip)
        code_llr = np.clip(dec.code_llr, -llr_clip, llr_clip)
        dec_ext = np.clip(code_llr - det_ext_deint, -llr_clip, llr_clip)
        detector_prior = dec_ext[interleaver]

    hard_info = (info_llr < 0.0).astype(np.uint8)
    errors = int(np.sum(hard_info != info_bits))
    return errors, int(info_bits.size)


def _shortened_turbo_frame(info_bits, v, n0, nu, turbo_iters, llr_clip,
                           interleaver, rng):
    encoded = conv_encode_75(info_bits)
    inv_perm = np.argsort(interleaver)
    int_encoded = encoded[interleaver]
    symbols = bpsk_modulate(int_encoded)
    # Paper [26] Forney model: y = V*x + white_noise
    y = forney_channel(symbols, v, n0=n0, rng=rng)

    n = len(y)
    n_fft = 1
    while n_fft < n + len(v):
        n_fft *= 2
    Z_w, g_diag, g_off, n_fft = compute_shortened_params_from_v(v, n0, nu, n_fft)
    z = apply_shortening_filter_fft(y, Z_w, n_fft)

    detector_prior = np.zeros(z.size, dtype=float)

    for _ in range(turbo_iters + 1):
        det = shortened_bcjr_bpsk(z, g_diag, g_off, nu=nu, la=detector_prior)
        det_llr = np.clip(det.llr, -llr_clip, llr_clip)
        det_ext = np.clip(det_llr - detector_prior, -llr_clip, llr_clip)
        det_ext_deint = det_ext[inv_perm]

        dec = conv_bcjr_decode(det_ext_deint)
        info_llr = np.clip(dec.info_llr, -llr_clip, llr_clip)
        code_llr = np.clip(dec.code_llr, -llr_clip, llr_clip)
        dec_ext = np.clip(code_llr - det_ext_deint, -llr_clip, llr_clip)
        detector_prior = dec_ext[interleaver]

    hard_info = (info_llr < 0.0).astype(np.uint8)
    errors = int(np.sum(hard_info != info_bits))
    return errors, int(info_bits.size)


# ---------------------------------------------------------------------------
# Curve runner
# ---------------------------------------------------------------------------

def _run_curve(label, snr_db_values, n_info, frames_per_snr, frame_fn, rng,
               min_errors=100, max_frames=500):
    rows = []
    for snr_db in snr_db_values:
        t0 = time.time()
        bit_errors = 0
        total_bits = 0
        for frame_idx in range(frames_per_snr):
            info_bits = rng.integers(0, 2, size=n_info, dtype=np.uint8)
            errs, total = frame_fn(info_bits, snr_db, rng)
            bit_errors += errs
            total_bits += total
            if bit_errors >= min_errors and frame_idx >= 2:
                break
            if frame_idx + 1 >= max_frames:
                break
        ber = bit_errors / total_bits if total_bits > 0 else 0.0
        elapsed = time.time() - t0
        rows.append({"snr_db": snr_db, "bit_errors": bit_errors,
                      "total_bits": total_bits, "ber": ber})
        print(f"  [{label}] Eb/N0={snr_db:5.2f} dB  BER={ber:.6g}  "
              f"({bit_errors}/{total_bits})  {elapsed:.1f}s")
    return rows


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Turbo equalization BER comparison")
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--tau", type=float, default=0.5, choices=[0.5, 0.35])
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--seed", type=int, default=20260529)
    args = parser.parse_args()

    tau = args.tau
    rolloff = 0.3
    pulse_span = 15
    sps = 128
    code_rate = 0.5
    turbo_iters = 5
    llr_clip = 20.0

    if tau == 0.5:
        isi_len = 7
        snr_range = list(np.arange(2.0, 6.1, 0.5))
    else:
        isi_len = 10
        snr_range = list(np.arange(2.5, 6.6, 0.5))

    if args.quick:
        K = 500
        frames_per_snr = 5
        min_errors = 10
        max_frames = 10
        snr_range = snr_range[::2]
    else:
        K = 6000
        frames_per_snr = 100
        min_errors = 200
        max_frames = 500

    output_dir = args.output_dir or f"results/comparison/turbo_tau{int(tau*1000):04d}"
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    t, h = generate_rrc(beta=rolloff, span=pulse_span, sps=sps)
    g = compute_g(t, h, tau=tau, isi_len=isi_len)
    v = min_phase_from_pulse(t, h, tau=tau, g0=g[0])

    interleaver = _make_interleaver(2 * K, seed=args.seed + 1000)
    rng = np.random.default_rng(args.seed)

    es_n0_offset = 10.0 * np.log10(code_rate)

    g_str = ", ".join(f"g[{k}]={val:.4f}" for k, val in sorted(g.items()) if abs(val) > 1e-6)
    print(f"Turbo equalization comparison: tau={tau}, isi_len={isi_len}, K={K}")
    print(f"g taps: {g_str}")
    print()

    all_results = {}

    # No-ISI baseline
    print("=== No-ISI baseline ===")
    all_results["no_isi"] = _run_curve(
        "no-ISI", snr_range, K, frames_per_snr,
        lambda bits, eb_n0, r: _noisi_frame(bits, eb_n0 + es_n0_offset, r),
        rng, min_errors=min_errors, max_frames=max_frames,
    )

    # Ungerboeck M-BCJR
    for M in [4, 8]:
        label = f"ung_M{M}"
        print(f"\n=== Ungerboeck M-BCJR M={M} ===")
        all_results[label] = _run_curve(
            label, snr_range, K, frames_per_snr,
            lambda bits, eb_n0, r, m=M: _ungerboeck_turbo_frame(
                bits, g, 1.0 / (10.0 ** ((eb_n0 + es_n0_offset) / 10.0)),
                isi_len, m, future_len=3, turbo_iters=turbo_iters,
                llr_clip=llr_clip, interleaver=interleaver, rng=r),
            rng, min_errors=min_errors, max_frames=max_frames,
        )

    # Paper [14] Smoothed Backup M-BCJR
    for M in [4, 8]:
        label = f"prlja_M{M}"
        print(f"\n=== Paper [14] Backup M-BCJR M={M} ===")
        all_results[label] = _run_curve(
            label, snr_range, K, frames_per_snr,
            lambda bits, eb_n0, r, m=M: _prlja_turbo_frame(
                bits, v, 1.0 / (10.0 ** ((eb_n0 + es_n0_offset) / 10.0)),
                m, turbo_iters=turbo_iters, llr_clip=llr_clip,
                interleaver=interleaver, rng=r),
            rng, min_errors=min_errors, max_frames=max_frames,
        )

    # Paper [26] Shortened BCJR (Forney model)
    for nu in [2, 3]:
        if nu >= isi_len:
            continue
        label = f"cs_nu{nu}"
        print(f"\n=== Paper [26] Shortened BCJR nu={nu} ===")
        all_results[label] = _run_curve(
            label, snr_range, K, frames_per_snr,
            lambda bits, eb_n0, r, n=nu: _shortened_turbo_frame(
                bits, v, 1.0 / (10.0 ** ((eb_n0 + es_n0_offset) / 10.0)),
                n, turbo_iters=turbo_iters, llr_clip=llr_clip,
                interleaver=interleaver, rng=r),
            rng, min_errors=min_errors, max_frames=max_frames,
        )

    # ---- Save CSV ----
    csv_path = output_path / "turbo_ber.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["snr_db"] + list(all_results.keys()))
        for i, snr in enumerate(snr_range):
            row = [snr] + [f"{all_results[k][i]['ber']:.8g}" for k in all_results]
            writer.writerow(row)
    print(f"\nSaved CSV to {csv_path}")

    # ---- Save config ----
    config = {
        "tau": tau, "rolloff": rolloff, "pulse_span": pulse_span, "sps": sps,
        "isi_len": isi_len, "K": K, "frames_per_snr": frames_per_snr,
        "seed": args.seed, "snr_range": snr_range,
        "turbo_iters": turbo_iters, "llr_clip": llr_clip,
        "g": {str(k): val for k, val in sorted(g.items())},
        "v": v.tolist(),
    }
    with (output_path / "config.json").open("w", encoding="utf-8") as f:
        json.dump(config, f, indent=2)

    # ---- Save NPZ ----
    np.savez(
        output_path / "turbo_data.npz",
        snr_db=np.array(snr_range),
        **{k: np.array([r["ber"] for r in v]) for k, v in all_results.items()},
    )

    print(f"Results saved to {output_path}")
    print("Done.")


if __name__ == "__main__":
    main()
