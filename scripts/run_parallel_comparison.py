"""Highly parallel comparison experiment runner.

Uses multiprocessing to parallelize across all CPU cores.
Each (algorithm, SNR point) pair is a separate task.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

PROJECT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT))

from ftn.baselines.forney_model import forney_channel, min_phase_from_pulse
from ftn.baselines.prlja_mbcjr import prlja_mbcjr_bpsk, prlja_backup_mbcjr_bpsk
from ftn.baselines.shortened_bcjr import cs_equalizer_bpsk, shortened_bcjr_bpsk
from ftn.baselines.channel_shortening import (
    compute_shortened_params_from_v,
    apply_shortening_filter_fft,
)
from ftn.channel import ftn_awgn_channel
from ftn.coding.conv import conv_bcjr_decode, conv_encode_75
from ftn.equalizers.full_bcjr import full_bcjr_bpsk
from ftn.equalizers.ungerboeck_mbcjr import ungerboeck_mbcjr_bpsk
from ftn.modulation import bpsk_modulate
from ftn.pulse import compute_g, generate_rrc


# ---------------------------------------------------------------------------
# Shared channel setup (computed once, passed to workers)
# ---------------------------------------------------------------------------

def _setup(tau, isi_len):
    t, h = generate_rrc(beta=0.3, span=15, sps=128)
    g = compute_g(t, h, tau=tau, isi_len=isi_len)
    v = min_phase_from_pulse(t, h, tau=tau, g0=g[0])
    return g, v


# ---------------------------------------------------------------------------
# Worker: run one (algorithm, snr) point
# ---------------------------------------------------------------------------

def _run_point(args):
    """Run frames for one algorithm at one SNR point. Returns dict."""
    (algo_id, snr_db, n_bits, frames_per_snr, min_errors, max_frames,
     seed, g, v, isi_len, code_rate, turbo_iters, llr_clip) = args

    n0 = 1.0 / (10.0 ** (snr_db / 10.0))
    es_n0_offset = 10.0 * np.log10(code_rate) if code_rate < 1 else 0.0
    n0_es = 1.0 / (10.0 ** ((snr_db + es_n0_offset) / 10.0))
    rng = np.random.default_rng(seed)

    bit_errors = 0
    total_bits = 0

    if algo_id == "full_bcjr":
        for _ in range(frames_per_snr):
            x = rng.choice([-1.0, 1.0], size=n_bits)
            y = ftn_awgn_channel(x, g, n0=n0, rng=rng)
            result = full_bcjr_bpsk(y, g, n0=n0, isi_len=isi_len)
            bit_errors += int(np.sum(np.sign(result.llr) != x))
            total_bits += n_bits
            if bit_errors >= min_errors:
                break

    elif algo_id.startswith("ung_M"):
        M = int(algo_id.split("M")[1])
        initial_state = tuple(1.0 for _ in range(isi_len))
        for _ in range(frames_per_snr):
            x = rng.choice([-1.0, 1.0], size=n_bits)
            y = ftn_awgn_channel(x, g, n0=n0, rng=rng)
            result = ungerboeck_mbcjr_bpsk(
                y, g, n0=n0, isi_len=isi_len, m_states=M,
                future_len=3, initial_state=initial_state)
            bit_errors += int(np.sum(np.sign(result.llr) != x))
            total_bits += n_bits
            if bit_errors >= min_errors:
                break

    elif algo_id.startswith("prlja_M"):
        M = int(algo_id.split("M")[1])
        for _ in range(frames_per_snr):
            x = rng.choice([-1.0, 1.0], size=n_bits)
            y = forney_channel(x, v, n0=n0, rng=rng)
            result = prlja_mbcjr_bpsk(y, v, n0=n0, M=M)
            bit_errors += int(np.sum(np.sign(result.llr) != x))
            total_bits += n_bits
            if bit_errors >= min_errors:
                break

    elif algo_id.startswith("cs_nu"):
        nu = int(algo_id.split("nu")[1])
        n_sig = n_bits
        n_fft = 1
        while n_fft < n_sig + len(v):
            n_fft *= 2
        Z_w, g_diag, g_off, _ = compute_shortened_params_from_v(v, n0, nu, n_fft)
        for _ in range(frames_per_snr):
            x = rng.choice([-1.0, 1.0], size=n_sig)
            y = forney_channel(x, v, n0=n0, rng=rng)
            z = apply_shortening_filter_fft(y, Z_w, n_fft)
            result = shortened_bcjr_bpsk(z, g_diag, g_off, nu=nu)
            bit_errors += int(np.sum(np.sign(result.llr) != x))
            total_bits += n_sig
            if bit_errors >= min_errors:
                break

    elif algo_id == "no_isi":
        # Turbo no-ISI baseline
        for _ in range(frames_per_snr):
            info_bits = rng.integers(0, 2, size=n_bits, dtype=np.uint8)
            encoded = conv_encode_75(info_bits)
            symbols = bpsk_modulate(encoded)
            noise = rng.standard_normal(symbols.size) * np.sqrt(n0_es / 2.0)
            y = symbols + noise
            code_llr = 4.0 * y / n0_es
            dec = conv_bcjr_decode(code_llr)
            hard_info = (dec.info_llr < 0.0).astype(np.uint8)
            bit_errors += int(np.sum(hard_info != info_bits))
            total_bits += int(info_bits.size)
            if bit_errors >= min_errors:
                break

    elif algo_id.startswith("turbo_ung_M"):
        M = int(algo_id.split("M")[1])
        interleaver = np.random.default_rng(seed + 1000).permutation(2 * n_bits)
        inv_perm = np.argsort(interleaver)
        initial_state = tuple(1.0 for _ in range(isi_len))
        for _ in range(frames_per_snr):
            info_bits = rng.integers(0, 2, size=n_bits, dtype=np.uint8)
            encoded = conv_encode_75(info_bits)
            int_encoded = encoded[interleaver]
            symbols = bpsk_modulate(int_encoded)
            y = ftn_awgn_channel(symbols, g, n0=n0_es, rng=rng)
            detector_prior = np.zeros(y.size, dtype=float)
            for _it in range(turbo_iters + 1):
                det = ungerboeck_mbcjr_bpsk(
                    y, g, n0=n0_es, la=detector_prior, isi_len=isi_len,
                    m_states=M, future_len=3, initial_state=initial_state)
                det_llr = np.clip(det.llr, -llr_clip, llr_clip)
                det_ext = np.clip(det_llr - detector_prior, -llr_clip, llr_clip)
                det_ext_deint = det_ext[inv_perm]
                dec = conv_bcjr_decode(det_ext_deint)
                info_llr = np.clip(dec.info_llr, -llr_clip, llr_clip)
                code_llr = np.clip(dec.code_llr, -llr_clip, llr_clip)
                dec_ext = np.clip(code_llr - det_ext_deint, -llr_clip, llr_clip)
                detector_prior = dec_ext[interleaver]
            hard_info = (info_llr < 0.0).astype(np.uint8)
            bit_errors += int(np.sum(hard_info != info_bits))
            total_bits += int(info_bits.size)
            if bit_errors >= min_errors:
                break

    elif algo_id.startswith("turbo_prlja_M"):
        M = int(algo_id.split("M")[1])
        interleaver = np.random.default_rng(seed + 1000).permutation(2 * n_bits)
        inv_perm = np.argsort(interleaver)
        for _ in range(frames_per_snr):
            info_bits = rng.integers(0, 2, size=n_bits, dtype=np.uint8)
            encoded = conv_encode_75(info_bits)
            int_encoded = encoded[interleaver]
            symbols = bpsk_modulate(int_encoded)
            y = forney_channel(symbols, v, n0=n0_es, rng=rng)
            detector_prior = np.zeros(y.size, dtype=float)
            for _it in range(turbo_iters + 1):
                det = prlja_backup_mbcjr_bpsk(
                    y, v, n0=n0_es, M=M, M_B=2, smooth=True, la=detector_prior)
                det_llr = np.clip(det.llr, -llr_clip, llr_clip)
                det_ext = np.clip(det_llr - detector_prior, -llr_clip, llr_clip)
                det_ext_deint = det_ext[inv_perm]
                dec = conv_bcjr_decode(det_ext_deint)
                info_llr = np.clip(dec.info_llr, -llr_clip, llr_clip)
                code_llr = np.clip(dec.code_llr, -llr_clip, llr_clip)
                dec_ext = np.clip(code_llr - det_ext_deint, -llr_clip, llr_clip)
                detector_prior = dec_ext[interleaver]
            hard_info = (info_llr < 0.0).astype(np.uint8)
            bit_errors += int(np.sum(hard_info != info_bits))
            total_bits += int(info_bits.size)
            if bit_errors >= min_errors:
                break

    elif algo_id.startswith("turbo_cs_nu"):
        nu = int(algo_id.split("nu")[1])
        interleaver = np.random.default_rng(seed + 1000).permutation(2 * n_bits)
        inv_perm = np.argsort(interleaver)
        n_fft = 1
        while n_fft < 2 * n_bits + 2 * len(v):
            n_fft *= 2
        Z_w, g_diag, g_off, _ = compute_shortened_params_from_v(v, n0_es, nu, n_fft)
        for _ in range(frames_per_snr):
            info_bits = rng.integers(0, 2, size=n_bits, dtype=np.uint8)
            encoded = conv_encode_75(info_bits)
            int_encoded = encoded[interleaver]
            symbols = bpsk_modulate(int_encoded)
            y = forney_channel(symbols, v, n0=n0_es, rng=rng)
            z = apply_shortening_filter_fft(y, Z_w, n_fft)
            detector_prior = np.zeros(z.size, dtype=float)
            for _it in range(turbo_iters + 1):
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
            bit_errors += int(np.sum(hard_info != info_bits))
            total_bits += int(info_bits.size)
            if bit_errors >= min_errors:
                break

    ber = bit_errors / total_bits if total_bits > 0 else 0.0
    return {"algo": algo_id, "snr_db": snr_db, "bit_errors": bit_errors,
            "total_bits": total_bits, "ber": ber}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--mode", default="all", choices=["uncoded", "turbo", "all"])
    parser.add_argument("--tau", nargs="*", default=[0.5, 0.35], type=float)
    parser.add_argument("--seed", type=int, default=20260529)
    parser.add_argument("--workers", type=int, default=None)
    args = parser.parse_args()

    n_workers = args.workers or min(os.cpu_count() or 4, 80)
    print(f"Using {n_workers} workers")

    for tau in args.tau:
        isi_len = 7 if tau == 0.5 else 10
        g, v = _setup(tau, isi_len)

        if args.quick:
            n_bits, frames, min_err, max_frames = 500, 5, 10, 10
        else:
            n_bits, frames, min_err, max_frames = 6000, 100, 200, 500

        if tau == 0.5:
            snr_range = list(np.arange(2.0, 6.1, 0.5))
        else:
            snr_range = list(np.arange(2.5, 6.6, 0.5))

        if args.quick:
            snr_range = snr_range[::2]

        tasks = []

        if args.mode in ("uncoded", "all"):
            uncoded_algos = ["full_bcjr", "ung_M4", "ung_M8", "prlja_M4", "prlja_M8"]
            for nu in [2, 3, 5]:
                if nu < isi_len:
                    uncoded_algos.append(f"cs_nu{nu}")
            for algo in uncoded_algos:
                for snr in snr_range:
                    tasks.append((algo, snr, n_bits, frames, min_err, max_frames,
                                  args.seed, g, v, isi_len, 1.0, 5, 20.0))

        if args.mode in ("turbo", "all"):
            turbo_algos = ["no_isi", "turbo_ung_M4", "turbo_ung_M8",
                           "turbo_prlja_M4", "turbo_prlja_M8"]
            for nu in [2, 3]:
                if nu < isi_len:
                    turbo_algos.append(f"turbo_cs_nu{nu}")
            for algo in turbo_algos:
                for snr in snr_range:
                    tasks.append((algo, snr, n_bits, frames, min_err, max_frames,
                                  args.seed, g, v, isi_len, 0.5, 5, 20.0))

        print(f"\ntau={tau}: {len(tasks)} tasks, {len(snr_range)} SNR points")

        results = {}
        t0 = time.time()
        with ProcessPoolExecutor(max_workers=n_workers) as pool:
            futures = {pool.submit(_run_point, t): t for t in tasks}
            done_count = 0
            for future in as_completed(futures):
                res = future.result()
                algo = res["algo"]
                if algo not in results:
                    results[algo] = []
                results[algo].append(res)
                done_count += 1
                if done_count % 10 == 0 or done_count == len(tasks):
                    elapsed = time.time() - t0
                    print(f"  [{tau}] {done_count}/{len(tasks)} done "
                          f"({elapsed:.0f}s, {done_count/elapsed:.1f}/s)")

        # Sort results by SNR
        for algo in results:
            results[algo].sort(key=lambda r: r["snr_db"])

        # Save outputs
        for mode, algos in [("uncoded", [a for a in results if not a.startswith("turbo") and a != "no_isi"]),
                            ("turbo", [a for a in results if a.startswith("turbo") or a == "no_isi"])]:
            if not algos:
                continue
            out_dir = Path(f"results/comparison/{mode}_tau{int(tau*1000):04d}")
            out_dir.mkdir(parents=True, exist_ok=True)

            csv_path = out_dir / f"{mode}_ber.csv"
            with csv_path.open("w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["snr_db"] + algos)
                for i, snr in enumerate(snr_range):
                    row = [snr] + [f"{results[a][i]['ber']:.8g}" for a in algos]
                    writer.writerow(row)
            print(f"Saved {csv_path}")

            config = {
                "tau": tau, "isi_len": isi_len, "K": n_bits,
                "seed": args.seed, "snr_range": snr_range,
            }
            with (out_dir / "config.json").open("w") as f:
                json.dump(config, f, indent=2)

            np.savez(
                out_dir / f"{mode}_data.npz",
                snr_db=np.array(snr_range),
                **{a: np.array([r["ber"] for r in results[a]]) for a in algos},
            )

    # Generate plots (best-effort, may fail if matplotlib has issues)
    for tau in args.tau:
        for mode in ["uncoded", "turbo"]:
            cmd = [sys.executable, str(PROJECT / "scripts/plot_comparison.py"),
                   mode, "--tau", str(tau)]
            try:
                subprocess.run(cmd, cwd=str(PROJECT), timeout=30)
            except Exception as e:
                print(f"Plotting failed for {mode} tau={tau}: {e}")

    print(f"\nTotal time: {time.time() - t0:.0f}s")
    print("Done.")


if __name__ == "__main__":
    import subprocess
    main()
