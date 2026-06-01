"""Reproduce Fig. 8: BPSK, tau=0.35, (7,5) conv code, turbo equalization.

Curves:
  - no-ISI baseline: BPSK + (7,5) code in AWGN (tau=1, no ISI)
  - FTN M=8, L=5: proposed M-BCJR with 8 survivors, future_len=5
  - FTN M=8, L=7: proposed M-BCJR with 8 survivors, future_len=7

Expected: M=8,L=7 approaches no-ISI around Eb/N0=4.5 dB.

SNR convention: the paper's x-axis is Eb/N0. For rate-1/2 coding with
unit-energy BPSK symbols, Es/N0(dB) = Eb/N0(dB) - 3.01.
The script takes Eb/N0 as input and converts internally.
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
from ftn.coding.conv import conv_bcjr_decode, conv_encode_75
from ftn.equalizers.ungerboeck_mbcjr import ungerboeck_mbcjr_bpsk
from ftn.modulation import bpsk_modulate
from ftn.pulse import compute_g, generate_rrc


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_interleaver(n: int, seed: int) -> np.ndarray:
    return np.random.default_rng(seed).permutation(n)


def _noisi_frame(info_bits: np.ndarray, es_n0_db: float, rng: np.random.Generator) -> tuple[int, int]:
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


def _ftn_turbo_frame(
    info_bits: np.ndarray,
    g: dict[int, float],
    n0: float,
    isi_len: int,
    m_states: int,
    future_len: int,
    turbo_iters: int,
    llr_clip: float,
    interleaver: np.ndarray,
    rng: np.random.Generator,
) -> tuple[int, int]:
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


def _run_curve(
    label: str,
    snr_db_values: list[float],
    n_info: int,
    frames_per_snr: int,
    frame_fn,
    rng: np.random.Generator,
    min_errors: int = 100,
    max_frames: int = 500,
) -> list[dict]:
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
        rows.append({
            "snr_db": snr_db, "bit_errors": bit_errors,
            "total_bits": total_bits, "ber": ber,
        })
        print(f"  [{label}] SNR={snr_db:5.2f} dB  BER={ber:.6g}  "
              f"({bit_errors}/{total_bits})  {elapsed:.1f}s")
    return rows


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="Reproduce Fig. 8: BPSK tau=0.35")
    parser.add_argument("--quick", action="store_true",
                        help="Quick mode: K=500, few frames, coarse SNR grid")
    parser.add_argument("--output-dir", default="results/fig8_tau_035")
    parser.add_argument("--seed", type=int, default=20260508)
    parser.add_argument("--K", type=int, default=None, help="Override info bits per frame")
    parser.add_argument("--frames", type=int, default=None, help="Override max frames per SNR")
    args = parser.parse_args()

    if args.quick:
        K = 500
        frames_per_snr = 5
        snr_range = [2.5, 3.5, 4.5, 5.5, 6.5]
        min_errors = 10
    else:
        K = 6000
        frames_per_snr = 100
        snr_range = list(np.arange(2.5, 6.5, 0.25))
        min_errors = 200

    if args.K is not None:
        K = args.K
    if args.frames is not None:
        frames_per_snr = args.frames

    tau = 0.35
    rolloff = 0.3
    pulse_span = 15
    sps = 128
    isi_len = 10
    turbo_iters = 15
    llr_clip = 20.0
    code_rate = 0.5

    configs = [
        {"label": "M8_L5", "m_states": 8, "future_len": 5},
        {"label": "M8_L7", "m_states": 8, "future_len": 7},
    ]

    output_path = Path(args.output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    t, h = generate_rrc(beta=rolloff, span=pulse_span, sps=sps)
    g = compute_g(t, h, tau=tau, isi_len=isi_len)

    interleaver = _make_interleaver(2 * K, seed=args.seed + 1000)
    rng = np.random.default_rng(args.seed)

    # snr_range is Eb/N0 from the paper; convert to Es/N0 for the channel
    eb_n0_range = snr_range
    es_n0_offset = 10.0 * np.log10(code_rate)  # -3.01 dB

    print(f"Fig. 8 reproduction: tau={tau}, K={K}, isi_len={isi_len}, "
          f"turbo_iters={turbo_iters}")
    print(f"SNR range: Eb/N0 = {eb_n0_range[0]:.1f} .. {eb_n0_range[-1]:.1f} dB "
          f"(Es/N0 = {eb_n0_range[0]+es_n0_offset:.1f} .. {eb_n0_range[-1]+es_n0_offset:.1f} dB)")
    g_str = ", ".join(f"g[{k}]={v:.4f}" for k, v in sorted(g.items()) if abs(v) > 1e-6)
    print(f"g taps: {g_str}")
    print()

    all_results: dict[str, list[dict]] = {}

    # No-ISI baseline (Eb/N0 → Es/N0 conversion)
    print("=== No-ISI baseline ===")
    all_results["no_isi"] = _run_curve(
        "no-ISI", eb_n0_range, K, frames_per_snr,
        lambda bits, eb_n0, r: _noisi_frame(bits, eb_n0 + es_n0_offset, r),
        rng, min_errors=min_errors,
    )

    # FTN curves
    for cfg in configs:
        print(f"\n=== FTN {cfg['label']} ===")
        all_results[cfg["label"]] = _run_curve(
            cfg["label"], eb_n0_range, K, frames_per_snr,
            lambda bits, eb_n0, r, c=cfg: _ftn_turbo_frame(
                bits, g, 1.0 / (10.0 ** ((eb_n0 + es_n0_offset) / 10.0)),
                isi_len, c["m_states"], c["future_len"],
                turbo_iters, llr_clip, interleaver, r,
            ),
            rng, min_errors=min_errors,
        )

    # ---- Save CSV ----
    csv_path = output_path / "fig8_results.csv"
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
        "isi_len": isi_len, "turbo_iters": turbo_iters, "llr_clip": llr_clip,
        "K": K, "frames_per_snr": frames_per_snr, "seed": args.seed,
        "snr_range": snr_range,
        "g": {str(k): v for k, v in sorted(g.items())},
        "configs": configs,
    }
    with (output_path / "config.json").open("w", encoding="utf-8") as f:
        json.dump(config, f, indent=2)

    # ---- Save NPZ ----
    np.savez(
        output_path / "fig8_data.npz",
        snr_db=np.array(snr_range),
        **{k: np.array([r["ber"] for r in v]) for k, v in all_results.items()},
    )

    print(f"Results saved to {output_path}")
    print("Done.")


if __name__ == "__main__":
    main()
