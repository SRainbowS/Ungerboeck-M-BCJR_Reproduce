"""Reproduce Fig. 10: BPSK, turbo code (R=1/3), turbo equalization.

Curves:
  - no-ISI baseline: BPSK + turbo code in AWGN
  - FTN tau=1.0 (no ISI), M=8, L=5
  - FTN tau=2/3, M=8, L=5
  - FTN tau=0.5, M=8, L=5

SNR conversion: Es/N0 = Eb/N0 + 10*log10(R) = Eb/N0 - 4.77 dB for R=1/3.
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

from ftn.channel import ftn_awgn_channel
from ftn.coding.turbo import (
    s_random_interleaver,
    turbo_decode,
    turbo_encode,
)
from ftn.equalizers.ungerboeck_mbcjr import ungerboeck_mbcjr_bpsk
from ftn.modulation import bpsk_modulate
from ftn.pulse import compute_g, generate_rrc


def _noisi_frame_turbo(info_bits, interleaver, es_n0_db, rng):
    n0 = 1.0 / (10.0 ** (es_n0_db / 10.0))
    enc = turbo_encode(info_bits, interleaver)
    symbols = 1.0 - 2.0 * enc.coded_bits.astype(float)
    noise = rng.standard_normal(symbols.size) * np.sqrt(n0 / 2.0)
    y = symbols + noise
    llr = 4.0 * y / n0
    n = info_bits.size + 5  # K + mu1 + mu2
    sys_llr = llr[:n]
    p1_llr = llr[n:2 * n]
    p2_llr = llr[2 * n:3 * n]
    result = turbo_decode(sys_llr, p1_llr, p2_llr, interleaver, max_iterations=50)
    hard = (result.info_llr < 0).astype(np.uint8)
    errors = int(np.sum(hard != info_bits))
    return errors, int(info_bits.size)


def _ftn_turbo_eq_frame(
    info_bits, g, n0, isi_len, m_states, future_len,
    turbo_iters, llr_clip, interleaver, rng,
):
    enc = turbo_encode(info_bits, interleaver)
    K = info_bits.size
    n_total = K + 5  # K + mu1 + mu2

    # ALL coded bits go through FTN channel
    symbols = 1.0 - 2.0 * enc.coded_bits.astype(float)
    y = ftn_awgn_channel(symbols, g, n0=n0, rng=rng)

    initial_state = tuple(1.0 for _ in range(isi_len))
    detector_prior = np.zeros(y.size, dtype=float)

    for it in range(turbo_iters):
        # Detector: equalize FTN ISI
        det = ungerboeck_mbcjr_bpsk(
            y, g, n0=n0, la=detector_prior,
            isi_len=isi_len, m_states=m_states,
            future_len=future_len, initial_state=initial_state,
        )
        det_llr = np.clip(det.llr, -llr_clip, llr_clip)
        det_ext = np.clip(det_llr - detector_prior, -llr_clip, llr_clip)

        # Demux detector extrinsic into sys/p1/p2
        sys_llr = det_ext[:n_total]
        p1_llr = det_ext[n_total:2 * n_total]
        p2_llr = det_ext[2 * n_total:3 * n_total]

        # Turbo decode (1 internal iteration per outer iteration)
        result = turbo_decode(sys_llr, p1_llr, p2_llr, interleaver, max_iterations=1)

        # Build detector prior from decoder's systematic extrinsic
        detector_prior = np.zeros(y.size, dtype=float)
        if result.sys_extrinsic is not None:
            detector_prior[:K] = np.clip(result.sys_extrinsic, -llr_clip, llr_clip)

    # Final full turbo decode
    det = ungerboeck_mbcjr_bpsk(
        y, g, n0=n0, la=detector_prior,
        isi_len=isi_len, m_states=m_states,
        future_len=future_len, initial_state=initial_state,
    )
    det_llr = np.clip(det.llr, -llr_clip, llr_clip)
    det_ext = np.clip(det_llr - detector_prior, -llr_clip, llr_clip)
    sys_llr = det_ext[:n_total]
    p1_llr = det_ext[n_total:2 * n_total]
    p2_llr = det_ext[2 * n_total:3 * n_total]
    result = turbo_decode(sys_llr, p1_llr, p2_llr, interleaver, max_iterations=50)
    hard = (result.info_llr < 0).astype(np.uint8)
    errors = int(np.sum(hard != info_bits))
    return errors, int(info_bits.size)


def _run_curve(label, snr_db_values, n_info, frames_per_snr,
               frame_fn, rng, min_errors=100, max_frames=500):
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
        print(f"  [{label}] SNR={snr_db:5.2f} dB  BER={ber:.6g}  "
              f"({bit_errors}/{total_bits})  {elapsed:.1f}s")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Reproduce Fig. 10: Turbo-coded BPSK")
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--output-dir", default="results/fig10")
    parser.add_argument("--seed", type=int, default=20260511)
    args = parser.parse_args()

    if args.quick:
        K = 500
        frames_per_snr = 3
        snr_range = [0.0, 1.0, 2.0, 3.0]
        min_errors = 10
    else:
        K = 21842
        frames_per_snr = 50
        snr_range = list(np.arange(-0.5, 3.6, 0.25))
        min_errors = 100

    rolloff = 0.3
    pulse_span = 15
    sps = 128
    isi_len = 3
    m_states = 8
    future_len = 5
    turbo_iters = 50
    llr_clip = 20.0
    code_rate = 1.0 / 3.0

    es_n0_offset = 10.0 * np.log10(code_rate)

    tau_list = [1.0, 2 / 3, 0.5]
    output_path = Path(args.output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    interleaver = s_random_interleaver(K, seed=args.seed + 100)

    print(f"Fig. 10 reproduction: K={K}, R=1/3, isi_len={isi_len}, M={m_states}, L={future_len}")
    print(f"SNR range: Eb/N0 = {snr_range[0]:.1f} .. {snr_range[-1]:.1f} dB")
    print()

    all_results = {}

    # No-ISI baseline
    print("=== No-ISI baseline ===")
    all_results["no_isi"] = _run_curve(
        "no-ISI", snr_range, K, frames_per_snr,
        lambda bits, eb_n0, r: _noisi_frame_turbo(
            bits, interleaver, eb_n0 + es_n0_offset, r),
        rng, min_errors=min_errors,
    )

    # FTN curves
    for tau in tau_list:
        t_h, h_h = generate_rrc(beta=rolloff, span=pulse_span, sps=sps)
        g = compute_g(t_h, h_h, tau=tau, isi_len=isi_len)
        label = f"tau={tau:.2f}"
        print(f"\n=== FTN {label} ===")
        all_results[f"ftn_tau{tau}"] = _run_curve(
            label, snr_range, K, frames_per_snr,
            lambda bits, eb_n0, r, _g=g: _ftn_turbo_eq_frame(
                bits, _g, 1.0 / (10.0 ** ((eb_n0 + es_n0_offset) / 10.0)),
                isi_len, m_states, future_len, turbo_iters, llr_clip,
                interleaver, r,
            ),
            rng, min_errors=min_errors,
        )

    # Save CSV
    csv_path = output_path / "fig10_results.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["snr_db"] + list(all_results.keys()))
        for i, snr in enumerate(snr_range):
            row = [snr] + [f"{all_results[k][i]['ber']:.8g}" for k in all_results]
            writer.writerow(row)

    config = {
        "tau_list": tau_list, "rolloff": rolloff, "pulse_span": pulse_span,
        "sps": sps, "isi_len": isi_len, "m_states": m_states,
        "future_len": future_len, "turbo_iters": turbo_iters,
        "K": K, "seed": args.seed, "snr_range": snr_range,
        "code_rate": code_rate,
    }
    with (output_path / "config.json").open("w", encoding="utf-8") as f:
        json.dump(config, f, indent=2)

    np.savez(
        output_path / "fig10_data.npz",
        snr_db=np.array(snr_range),
        **{k: np.array([r["ber"] for r in v]) for k, v in all_results.items()},
    )

    print(f"\nResults saved to {output_path}")
    print("Done.")


if __name__ == "__main__":
    main()
