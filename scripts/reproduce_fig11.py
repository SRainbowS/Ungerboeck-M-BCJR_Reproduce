"""Reproduce Fig. 11: 16-QAM, simplified M-BCJR + (7,5) conv code.

Curves:
  - no-ISI baseline: 16-QAM + (7,5) code in AWGN
  - FTN tau=2/3, M=4, L=3
  - FTN tau=0.8, M=4, L=3

SNR conversion: Es/N0 = Eb/N0 + 10*log10(R*log2(q)) = Eb/N0 + 3.01 dB
for R=1/2, q=16.
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

from ftn.channel import ftn_awgn_channel_complex
from ftn.coding.conv import conv_bcjr_decode, conv_encode_75
from ftn.equalizers.simplified_mbcjr import simplified_mbcjr_qam16
from ftn.modulation import bpsk_hard_bits_from_llr, qam16_modulate
from ftn.pulse import compute_g, generate_rrc


def _make_interleaver(n: int, seed: int) -> np.ndarray:
    return np.random.default_rng(seed).permutation(n)


def _noisi_frame_16qam(info_bits, interleaver, es_n0_db, rng):
    n0 = 1.0 / (10.0 ** (es_n0_db / 10.0))
    encoded = conv_encode_75(info_bits)
    n_code_bits = encoded.size
    if interleaver.size != n_code_bits:
        raise ValueError("interleaver length must match encoded codeword length.")
    inv_perm = np.argsort(interleaver)
    int_encoded = encoded[interleaver]
    pad = (4 - n_code_bits % 4) % 4
    if pad > 0:
        tx_bits = np.concatenate([int_encoded, np.zeros(pad, dtype=np.uint8)])
    else:
        tx_bits = int_encoded
    symbols = qam16_modulate(tx_bits)

    # AWGN: noise variance N0/2 per dimension for unit-energy constellation
    sigma2 = n0 / 2.0
    noise = rng.standard_normal(symbols.size) + 1j * rng.standard_normal(symbols.size)
    y = symbols + noise * np.sqrt(sigma2)

    # Soft demod: LLR per code bit using log-likelihood
    from ftn.modulation import qam16_constellation, qam16_demod_bit_llr
    symbols_ref, bit_map = qam16_constellation()
    n_sym = y.size
    log_post = np.zeros((n_sym, 16))
    for k in range(n_sym):
        for si in range(16):
            dist2 = abs(y[k] - symbols_ref[si]) ** 2
            log_post[k, si] = -dist2 / n0
    bit_llr = qam16_demod_bit_llr(log_post)

    code_llr = bit_llr[:n_code_bits][inv_perm]
    result = conv_bcjr_decode(code_llr)
    hard_info = (result.info_llr < 0.0).astype(np.uint8)
    errors = int(np.sum(hard_info != info_bits))
    return errors, int(info_bits.size)


def _ftn_turbo_frame_16qam(
    info_bits, g, n0, isi_len, m_states, future_len,
    turbo_iters, llr_clip, interleaver, rng,
):
    encoded = conv_encode_75(info_bits)
    n_code_bits = encoded.size
    if interleaver.size != n_code_bits:
        raise ValueError("interleaver length must match encoded codeword length.")
    inv_perm = np.argsort(interleaver)
    int_encoded = encoded[interleaver]
    # Pad to multiple of 4 for 16-QAM
    pad = (4 - n_code_bits % 4) % 4
    if pad > 0:
        tx_bits = np.concatenate([int_encoded, np.zeros(pad, dtype=np.uint8)])
    else:
        tx_bits = int_encoded
    n_sym = tx_bits.size // 4

    symbols = qam16_modulate(tx_bits)
    y = ftn_awgn_channel_complex(symbols, g, n0=n0, rng=rng)

    la = np.zeros(n_sym * 4, dtype=np.float64)

    for _ in range(turbo_iters + 1):
        det_result = simplified_mbcjr_qam16(
            y, g, n0=n0, la=la,
            isi_len=isi_len, m_states=m_states, future_len=future_len,
        )
        det_llr = np.clip(det_result.bit_llr, -llr_clip, llr_clip)
        det_ext = np.clip(det_llr - la, -llr_clip, llr_clip)
        code_llr_in = det_ext[:n_code_bits][inv_perm]

        dec = conv_bcjr_decode(code_llr_in)
        info_llr = np.clip(dec.info_llr, -llr_clip, llr_clip)
        code_llr = np.clip(dec.code_llr, -llr_clip, llr_clip)
        dec_ext = np.clip(code_llr - code_llr_in, -llr_clip, llr_clip)

        # Feed back as detector prior
        la = np.zeros(n_sym * 4, dtype=np.float64)
        la[:n_code_bits] = dec_ext[interleaver]

    hard_info = (info_llr < 0.0).astype(np.uint8)
    errors = int(np.sum(hard_info != info_bits))
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
    parser = argparse.ArgumentParser(description="Reproduce Fig. 11: 16-QAM simplified M-BCJR")
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--output-dir", default="results/fig11")
    parser.add_argument("--seed", type=int, default=20260509)
    args = parser.parse_args()

    if args.quick:
        K = 200
        frames_per_snr = 3
        snr_range = [8.0, 10.0, 12.0, 14.0]
        min_errors = 10
    else:
        K = 6000
        frames_per_snr = 100
        snr_range = list(np.arange(6.0, 16.1, 0.5))
        min_errors = 200

    rolloff = 0.3
    pulse_span = 15
    sps = 128
    isi_len = 3
    m_states = 4
    future_len = 3
    turbo_iters = 10
    llr_clip = 20.0
    code_rate = 0.5

    # Es/N0 = Eb/N0 + 10*log10(R*log2(16)) = Eb/N0 + 3.01 dB
    es_n0_offset = 10.0 * np.log10(code_rate * 4)

    tau_list = [2 / 3, 0.8]
    output_path = Path(args.output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    interleaver = _make_interleaver(2 * K, seed=args.seed + 1000)

    print(f"Fig. 11 reproduction: K={K}, isi_len={isi_len}, M={m_states}, L={future_len}")
    print(f"SNR range: Eb/N0 = {snr_range[0]:.1f} .. {snr_range[-1]:.1f} dB")
    print()

    all_results = {}

    # No-ISI baseline
    print("=== No-ISI baseline ===")
    all_results["no_isi"] = _run_curve(
        "no-ISI", snr_range, K, frames_per_snr,
        lambda bits, eb_n0, r: _noisi_frame_16qam(
            bits, interleaver, eb_n0 + es_n0_offset, r,
        ),
        rng, min_errors=min_errors,
    )

    # FTN curves
    for tau in tau_list:
        t_h, h_h = generate_rrc(beta=rolloff, span=pulse_span, sps=sps)
        g = compute_g(t_h, h_h, tau=tau, isi_len=isi_len)
        label = f"tau={tau:.4f}"
        print(f"\n=== FTN {label} ===")
        all_results[f"ftn_tau{tau}"] = _run_curve(
            label, snr_range, K, frames_per_snr,
            lambda bits, eb_n0, r, _g=g: _ftn_turbo_frame_16qam(
                bits, _g, 1.0 / (10.0 ** ((eb_n0 + es_n0_offset) / 10.0)),
                isi_len, m_states, future_len, turbo_iters, llr_clip,
                interleaver, r,
            ),
            rng, min_errors=min_errors,
        )

    # Save CSV
    csv_path = output_path / "fig11_results.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["snr_db"] + list(all_results.keys()))
        for i, snr in enumerate(snr_range):
            row = [snr] + [f"{all_results[k][i]['ber']:.8g}" for k in all_results]
            writer.writerow(row)

    # Save config
    config = {
        "tau_list": tau_list, "rolloff": rolloff, "pulse_span": pulse_span,
        "sps": sps, "isi_len": isi_len, "m_states": m_states,
        "future_len": future_len, "turbo_iters": turbo_iters,
        "K": K, "seed": args.seed, "snr_range": snr_range,
        "interleaver": "random permutation of encoded bits before 16-QAM mapper",
        "interleaver_seed": args.seed + 1000,
    }
    with (output_path / "config.json").open("w", encoding="utf-8") as f:
        json.dump(config, f, indent=2)

    # Save NPZ
    np.savez(
        output_path / "fig11_data.npz",
        snr_db=np.array(snr_range),
        **{k: np.array([r["ber"] for r in v]) for k, v in all_results.items()},
    )

    print(f"\nResults saved to {output_path}")
    print("Done.")


if __name__ == "__main__":
    main()
