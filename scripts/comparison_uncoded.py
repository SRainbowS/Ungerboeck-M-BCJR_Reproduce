"""Uncoded BPSK BER comparison: Ungerboeck M-BCJR vs Paper [14] vs Paper [26].

Produces BER vs Eb/N0 curves for all algorithms under identical FTN channel
conditions (same g-taps, same symbol sequences, same SNR).

Curves:
  - Full BCJR oracle (Ungerboeck model)
  - Ungerboeck M-BCJR, M=4 and M=8
  - Paper [14] Simple M-BCJR, M=4 and M=8
  - Paper [26] Shortened BCJR, nu=2, nu=3, nu=5

Output: results/comparison/uncoded_tau{05,035}_ber.csv + PNG
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

from ftn.baselines.forney_model import forney_channel, min_phase_from_pulse
from ftn.baselines.prlja_mbcjr import prlja_mbcjr_bpsk
from ftn.baselines.shortened_bcjr import cs_equalizer_bpsk
from ftn.channel import ftn_awgn_channel
from ftn.equalizers.full_bcjr import full_bcjr_bpsk
from ftn.equalizers.ungerboeck_mbcjr import ungerboeck_mbcjr_bpsk
from ftn.pulse import compute_g, generate_rrc


# ---------------------------------------------------------------------------
# Per-frame helpers for each equalizer
# ---------------------------------------------------------------------------

def _full_bcjr_frame(x, g, n0, rng):
    y = ftn_awgn_channel(x, g, n0=n0, rng=rng)
    isi_len = max(k for k in g if k > 0)
    result = full_bcjr_bpsk(y, g, n0=n0, isi_len=isi_len)
    hard = np.sign(result.llr)
    return int(np.sum(hard != x)), len(x)


def _ungerboeck_mbcjr_frame(x, g, n0, M, future_len, rng):
    y = ftn_awgn_channel(x, g, n0=n0, rng=rng)
    isi_len = max(k for k in g if k > 0)
    initial_state = tuple(1.0 for _ in range(isi_len))
    result = ungerboeck_mbcjr_bpsk(
        y, g, n0=n0, isi_len=isi_len, m_states=M,
        future_len=future_len, initial_state=initial_state,
    )
    hard = np.sign(result.llr)
    return int(np.sum(hard != x)), len(x)


def _prlja_mbcjr_frame(x, v, n0, M, rng):
    y = forney_channel(x, v, n0=n0, rng=rng)
    result = prlja_mbcjr_bpsk(y, v, n0=n0, M=M)
    hard = np.sign(result.llr)
    return int(np.sum(hard != x)), len(x)


def _shortened_bcjr_frame(x, v, n0, nu, rng):
    # Paper [26] Forney model: y = V*x + white_noise
    y = forney_channel(x, v, n0=n0, rng=rng)
    result = cs_equalizer_bpsk(y, v, n0=n0, nu=nu)
    hard = np.sign(result.llr)
    return int(np.sum(hard != x)), len(x)


# ---------------------------------------------------------------------------
# Curve runner
# ---------------------------------------------------------------------------

def _run_curve(label, snr_db_values, n_bits, frames_per_snr, frame_fn, rng,
               min_errors=100, max_frames=500):
    rows = []
    for snr_db in snr_db_values:
        t0 = time.time()
        bit_errors = 0
        total_bits = 0
        for frame_idx in range(frames_per_snr):
            x = rng.choice([-1.0, 1.0], size=n_bits)
            errs, total = frame_fn(x, snr_db, rng)
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
    parser = argparse.ArgumentParser(description="Uncoded BER comparison")
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--tau", type=float, default=0.5, choices=[0.5, 0.35])
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--seed", type=int, default=20260529)
    args = parser.parse_args()

    tau = args.tau
    rolloff = 0.3
    pulse_span = 15
    sps = 128

    if tau == 0.5:
        isi_len = 7
        snr_range = list(np.arange(2.0, 6.1, 0.5))
    else:
        isi_len = 10
        snr_range = list(np.arange(2.5, 6.6, 0.5))

    if args.quick:
        n_bits = 500
        frames_per_snr = 5
        min_errors = 10
        max_frames = 10
        snr_range = snr_range[::2]
    else:
        n_bits = 6000
        frames_per_snr = 100
        min_errors = 200
        max_frames = 500

    output_dir = args.output_dir or f"results/comparison/uncoded_tau{int(tau*1000):04d}"
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    t, h = generate_rrc(beta=rolloff, span=pulse_span, sps=sps)
    g = compute_g(t, h, tau=tau, isi_len=isi_len)
    v = min_phase_from_pulse(t, h, tau=tau, g0=g[0])

    rng = np.random.default_rng(args.seed)

    g_str = ", ".join(f"g[{k}]={val:.4f}" for k, val in sorted(g.items()) if abs(val) > 1e-6)
    print(f"Uncoded BER comparison: tau={tau}, isi_len={isi_len}, K={n_bits}")
    print(f"g taps: {g_str}")
    print(f"v (min-phase): {np.array2string(v, precision=4)}")
    print()

    all_results = {}

    # 1) Full BCJR oracle (Ungerboeck model)
    print("=== Full BCJR (Ungerboeck) ===")
    all_results["full_bcjr"] = _run_curve(
        "Full BCJR", snr_range, n_bits, frames_per_snr,
        lambda x, snr, r: _full_bcjr_frame(
            x, g, 1.0 / (10.0 ** (snr / 10.0)), r),
        rng, min_errors=min_errors, max_frames=max_frames,
    )

    # 2) Ungerboeck M-BCJR
    for M in [4, 8]:
        label = f"ung_M{M}"
        print(f"\n=== Ungerboeck M-BCJR M={M} ===")
        all_results[label] = _run_curve(
            label, snr_range, n_bits, frames_per_snr,
            lambda x, snr, r, m=M: _ungerboeck_mbcjr_frame(
                x, g, 1.0 / (10.0 ** (snr / 10.0)), m, future_len=3, rng=r),
            rng, min_errors=min_errors, max_frames=max_frames,
        )

    # 3) Paper [14] Simple M-BCJR (Forney model)
    for M in [4, 8]:
        label = f"prlja_M{M}"
        print(f"\n=== Paper [14] M-BCJR M={M} ===")
        all_results[label] = _run_curve(
            label, snr_range, n_bits, frames_per_snr,
            lambda x, snr, r, m=M: _prlja_mbcjr_frame(
                x, v, 1.0 / (10.0 ** (snr / 10.0)), m, r),
            rng, min_errors=min_errors, max_frames=max_frames,
        )

    # 4) Paper [26] Shortened BCJR (Forney model)
    for nu in [2, 3, 5]:
        if nu >= isi_len:
            continue  # nu >= L is equivalent to full BCJR, skip
        label = f"cs_nu{nu}"
        print(f"\n=== Paper [26] Shortened BCJR nu={nu} ===")
        all_results[label] = _run_curve(
            label, snr_range, n_bits, frames_per_snr,
            lambda x, snr, r, n=nu: _shortened_bcjr_frame(
                x, v, 1.0 / (10.0 ** (snr / 10.0)), n, r),
            rng, min_errors=min_errors, max_frames=max_frames,
        )

    # ---- Save CSV ----
    csv_path = output_path / "uncoded_ber.csv"
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
        "isi_len": isi_len, "K": n_bits, "frames_per_snr": frames_per_snr,
        "seed": args.seed, "snr_range": snr_range,
        "g": {str(k): val for k, val in sorted(g.items())},
        "v": v.tolist(),
    }
    with (output_path / "config.json").open("w", encoding="utf-8") as f:
        json.dump(config, f, indent=2)

    # ---- Save NPZ ----
    np.savez(
        output_path / "uncoded_data.npz",
        snr_db=np.array(snr_range),
        **{k: np.array([r["ber"] for r in v]) for k, v in all_results.items()},
    )

    print(f"Results saved to {output_path}")
    print("Done.")


if __name__ == "__main__":
    main()
