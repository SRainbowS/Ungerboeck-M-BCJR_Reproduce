"""Parallel Monte Carlo runner using multiprocessing — Windows/Linux portable.

Usage:
    python scripts/run_parallel.py fig7 --K 1000 --frames 8 --procs 8
    python scripts/run_parallel.py fig8 --K 500  --frames 5 --procs 4
    python scripts/run_parallel.py fig10 --K 2000 --frames 10 --procs 52
    python scripts/run_parallel.py fig11 --K 1000 --frames 10 --procs 52
    python scripts/run_parallel.py fig7  # full: K=6000, frames=100, procs=8
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

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ftn.channel import ftn_awgn_channel, ftn_awgn_channel_complex
from ftn.coding.conv import conv_bcjr_decode, conv_encode_75
from ftn.coding.turbo import (
    s_random_interleaver,
    turbo_decode,
    turbo_encode,
)
from ftn.equalizers.simplified_mbcjr import simplified_mbcjr_qam16
from ftn.equalizers.ungerboeck_mbcjr import ungerboeck_mbcjr_bpsk
from ftn.modulation import (
    bpsk_modulate,
    qam16_constellation,
    qam16_demod_bit_llr,
    qam16_modulate,
)
from ftn.pulse import compute_g, generate_rrc


# ── Figure parameters ───────────────────────────────────────────

FIGURE_PARAMS = {
    "fig7": {
        "tau": 0.5, "isi_len": 7, "turbo_iters": 5,
        "eb_n0_range": np.arange(2.0, 6.01, 0.25),
        "configs": [
            {"label": "M2_L3", "m_states": 2, "future_len": 3},
            {"label": "M2_L5", "m_states": 2, "future_len": 5},
        ],
        "output_dir": "results/fig7_tau_05",
    },
    "fig8": {
        "tau": 0.35, "isi_len": 10, "turbo_iters": 15,
        "eb_n0_range": np.arange(2.5, 6.51, 0.25),
        "configs": [
            {"label": "M8_L5", "m_states": 8, "future_len": 5},
            {"label": "M8_L7", "m_states": 8, "future_len": 7},
        ],
        "output_dir": "results/fig8_tau_035",
    },
    "fig10": {
        "tau_list": [1.0, 2 / 3, 0.5], "isi_len": 3,
        "m_states": 8, "future_len": 5, "turbo_iters": 50,
        "code_rate": 1.0 / 3.0,
        "eb_n0_range": np.arange(-0.5, 3.51, 0.25),
        "output_dir": "results/fig10",
    },
    "fig11": {
        "tau_list": [2 / 3, 0.8], "isi_len": 3,
        "m_states": 4, "future_len": 3, "turbo_iters": 10,
        "code_rate": 0.5,
        "eb_n0_range": np.arange(6.0, 16.01, 0.5),
        "output_dir": "results/fig11",
    },
}

CODE_RATE_CONV = 0.5
ES_N0_OFFSET_CONV = 10.0 * np.log10(CODE_RATE_CONV)


# ── Fig 7/8: BPSK + (7,5) conv code ────────────────────────────

def _run_one_task_fig78(args: tuple) -> dict:
    """Worker: run all frames for one (eb_n0, config) combination."""
    (eb_n0, label, m_states, future_len, tau, isi_len, turbo_iters,
     K, max_frames, min_errors, seed) = args

    es_n0 = eb_n0 + ES_N0_OFFSET_CONV
    n0 = 1.0 / (10.0 ** (es_n0 / 10.0))

    rng_local = np.random.default_rng(seed)
    t, h = generate_rrc(beta=0.3, span=15, sps=128)
    g = compute_g(t, h, tau=tau, isi_len=isi_len)
    interleaver = np.random.default_rng(seed + 1000).permutation(2 * K)
    inv_perm = np.argsort(interleaver)
    initial_state = tuple(1.0 for _ in range(isi_len))
    llr_clip = 20.0

    bit_errors = 0
    total_bits = 0
    for fi in range(max_frames):
        info_bits = rng_local.integers(0, 2, size=K, dtype=np.uint8)
        encoded = conv_encode_75(info_bits)
        int_encoded = encoded[interleaver]
        symbols = bpsk_modulate(int_encoded)
        y = ftn_awgn_channel(symbols, g, n0=n0, rng=rng_local)

        detector_prior = np.zeros(y.size, dtype=float)
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

        hard_info = (info_llr < 0).astype(np.uint8)
        bit_errors += int(np.sum(hard_info != info_bits))
        total_bits += int(info_bits.size)
        if bit_errors >= min_errors and fi >= 3:
            break

    ber = bit_errors / total_bits if total_bits > 0 else 0.0
    return {
        "eb_n0_db": eb_n0, "label": label,
        "m_states": m_states, "future_len": future_len,
        "bit_errors": bit_errors, "total_bits": total_bits, "ber": ber,
    }


def _run_noisi_baseline_fig78(
    eb_n0_range: list[float], K: int, max_frames: int,
    min_errors: int, seed: int,
) -> list[dict]:
    """No-ISI baseline for fig7/fig8 — fast, runs in main process."""
    rng = np.random.default_rng(seed)
    results = []
    for eb_n0 in eb_n0_range:
        es_n0 = eb_n0 + ES_N0_OFFSET_CONV
        n0 = 1.0 / (10.0 ** (es_n0 / 10.0))
        bit_errors = 0
        total_bits = 0
        for _ in range(max_frames):
            info = rng.integers(0, 2, size=K, dtype=np.uint8)
            enc = conv_encode_75(info)
            sym = bpsk_modulate(enc)
            y = sym + rng.standard_normal(sym.size) * np.sqrt(n0 / 2.0)
            code_llr = 4.0 * y / n0
            dec = conv_bcjr_decode(code_llr)
            hard = (dec.info_llr < 0).astype(np.uint8)
            bit_errors += int(np.sum(hard != info))
            total_bits += K
            if bit_errors >= min_errors:
                break
        ber = bit_errors / total_bits if total_bits > 0 else 0.0
        results.append({
            "eb_n0_db": eb_n0, "label": "no_isi",
            "bit_errors": bit_errors, "total_bits": total_bits, "ber": ber,
        })
        print(f"  [no-ISI] Eb/N0={eb_n0:.2f} dB  BER={ber:.6g}  "
              f"({bit_errors}/{total_bits})")
    return results


# ── Fig 10: BPSK + Turbo R=1/3 + turbo equalization ────────────

def _run_one_task_fig10(args: tuple) -> dict:
    """Worker: run all frames for one (eb_n0, tau) combination."""
    (eb_n0, label, tau, isi_len, m_states, future_len, turbo_iters,
     K, max_frames, min_errors, seed, code_rate) = args

    es_n0_offset = 10.0 * np.log10(code_rate)
    es_n0 = eb_n0 + es_n0_offset
    n0 = 1.0 / (10.0 ** (es_n0 / 10.0))

    rng_local = np.random.default_rng(seed)
    t, h = generate_rrc(beta=0.3, span=15, sps=128)
    g = compute_g(t, h, tau=tau, isi_len=isi_len)
    interleaver = s_random_interleaver(K, seed=seed + 100)

    llr_clip = 20.0
    initial_state = tuple(1.0 for _ in range(isi_len))
    n_total = K + 5  # K + mu1 + mu2

    bit_errors = 0
    total_bits = 0
    for fi in range(max_frames):
        info_bits = rng_local.integers(0, 2, size=K, dtype=np.uint8)
        enc = turbo_encode(info_bits, interleaver)
        symbols = 1.0 - 2.0 * enc.coded_bits.astype(float)
        y = ftn_awgn_channel(symbols, g, n0=n0, rng=rng_local)

        detector_prior = np.zeros(y.size, dtype=float)

        for _it in range(turbo_iters):
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

            result = turbo_decode(
                sys_llr, p1_llr, p2_llr, interleaver, max_iterations=1,
            )

            detector_prior = np.zeros(y.size, dtype=float)
            if result.sys_extrinsic is not None:
                detector_prior[:K] = np.clip(
                    result.sys_extrinsic, -llr_clip, llr_clip,
                )

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
        result = turbo_decode(
            sys_llr, p1_llr, p2_llr, interleaver, max_iterations=50,
        )

        hard = (result.info_llr < 0).astype(np.uint8)
        bit_errors += int(np.sum(hard != info_bits))
        total_bits += int(info_bits.size)
        if bit_errors >= min_errors and fi >= 3:
            break

    ber = bit_errors / total_bits if total_bits > 0 else 0.0
    return {
        "eb_n0_db": eb_n0, "label": label,
        "bit_errors": bit_errors, "total_bits": total_bits, "ber": ber,
    }


def _run_noisi_baseline_fig10(
    eb_n0_range: list[float], K: int, max_frames: int,
    min_errors: int, seed: int, code_rate: float, interleaver: np.ndarray,
) -> list[dict]:
    """No-ISI baseline for fig10: BPSK + turbo R=1/3 in AWGN."""
    es_n0_offset = 10.0 * np.log10(code_rate)
    rng = np.random.default_rng(seed)
    results = []
    for eb_n0 in eb_n0_range:
        es_n0 = eb_n0 + es_n0_offset
        n0 = 1.0 / (10.0 ** (es_n0 / 10.0))
        bit_errors = 0
        total_bits = 0
        for _ in range(max_frames):
            info = rng.integers(0, 2, size=K, dtype=np.uint8)
            enc = turbo_encode(info, interleaver)
            symbols = 1.0 - 2.0 * enc.coded_bits.astype(float)
            noise = rng.standard_normal(symbols.size) * np.sqrt(n0 / 2.0)
            y = symbols + noise
            llr = 4.0 * y / n0
            n_total = K + 5
            sys_llr = llr[:n_total]
            p1_llr = llr[n_total:2 * n_total]
            p2_llr = llr[2 * n_total:3 * n_total]
            result = turbo_decode(
                sys_llr, p1_llr, p2_llr, interleaver, max_iterations=50,
            )
            hard = (result.info_llr < 0).astype(np.uint8)
            bit_errors += int(np.sum(hard != info))
            total_bits += K
            if bit_errors >= min_errors:
                break
        ber = bit_errors / total_bits if total_bits > 0 else 0.0
        results.append({
            "eb_n0_db": eb_n0, "label": "no_isi",
            "bit_errors": bit_errors, "total_bits": total_bits, "ber": ber,
        })
        print(f"  [no-ISI] Eb/N0={eb_n0:.2f} dB  BER={ber:.6g}  "
              f"({bit_errors}/{total_bits})")
    return results


# ── Fig 11: 16-QAM + (7,5) conv code + simplified M-BCJR ───────

def _run_one_task_fig11(args: tuple) -> dict:
    """Worker: run all frames for one (eb_n0, tau) combination."""
    (eb_n0, label, tau, isi_len, m_states, future_len, turbo_iters,
     K, max_frames, min_errors, seed, code_rate) = args

    # Es/N0 = Eb/N0 + 10*log10(R*log2(q)) = Eb/N0 + 3.01 dB
    es_n0_offset = 10.0 * np.log10(code_rate * 4)
    n0 = 1.0 / (10.0 ** ((eb_n0 + es_n0_offset) / 10.0))

    rng_local = np.random.default_rng(seed)
    t, h = generate_rrc(beta=0.3, span=15, sps=128)
    g = compute_g(t, h, tau=tau, isi_len=isi_len)
    interleaver = np.random.default_rng(seed + 1000).permutation(2 * K)
    inv_perm = np.argsort(interleaver)
    llr_clip = 20.0

    bit_errors = 0
    total_bits = 0
    for fi in range(max_frames):
        info_bits = rng_local.integers(0, 2, size=K, dtype=np.uint8)
        encoded = conv_encode_75(info_bits)
        n_code_bits = encoded.size
        int_encoded = encoded[interleaver]
        pad = (4 - n_code_bits % 4) % 4
        if pad > 0:
            tx_bits = np.concatenate([int_encoded, np.zeros(pad, dtype=np.uint8)])
        else:
            tx_bits = int_encoded
        n_sym = tx_bits.size // 4

        symbols = qam16_modulate(tx_bits)
        y = ftn_awgn_channel_complex(symbols, g, n0=n0, rng=rng_local)

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

            la = np.zeros(n_sym * 4, dtype=np.float64)
            la[:n_code_bits] = dec_ext[interleaver]

        hard_info = (info_llr < 0.0).astype(np.uint8)
        bit_errors += int(np.sum(hard_info != info_bits))
        total_bits += int(info_bits.size)
        if bit_errors >= min_errors and fi >= 3:
            break

    ber = bit_errors / total_bits if total_bits > 0 else 0.0
    return {
        "eb_n0_db": eb_n0, "label": label,
        "bit_errors": bit_errors, "total_bits": total_bits, "ber": ber,
    }


def _run_noisi_baseline_fig11(
    eb_n0_range: list[float], K: int, max_frames: int,
    min_errors: int, seed: int, code_rate: float,
) -> list[dict]:
    """No-ISI baseline for fig11: 16-QAM + (7,5) code in AWGN."""
    es_n0_offset = 10.0 * np.log10(code_rate * 4)
    rng = np.random.default_rng(seed)
    symbols_ref, bit_map = qam16_constellation()
    interleaver = np.random.default_rng(seed + 1000).permutation(2 * K)
    inv_perm = np.argsort(interleaver)
    results = []
    for eb_n0 in eb_n0_range:
        es_n0 = eb_n0 + es_n0_offset
        n0 = 1.0 / (10.0 ** (es_n0 / 10.0))
        bit_errors = 0
        total_bits = 0
        for _ in range(max_frames):
            info = rng.integers(0, 2, size=K, dtype=np.uint8)
            enc = conv_encode_75(info)
            n_code_bits = enc.size
            int_enc = enc[interleaver]
            pad = (4 - n_code_bits % 4) % 4
            if pad > 0:
                tx_bits = np.concatenate([int_enc, np.zeros(pad, dtype=np.uint8)])
            else:
                tx_bits = int_enc
            symbols = qam16_modulate(tx_bits)

            sigma2 = n0 / 2.0
            noise = (rng.standard_normal(symbols.size)
                     + 1j * rng.standard_normal(symbols.size))
            y = symbols + noise * np.sqrt(sigma2)

            # Vectorized soft demod
            n_sym = y.size
            log_post = np.zeros((n_sym, 16))
            for si in range(16):
                log_post[:, si] = -np.abs(y - symbols_ref[si]) ** 2 / n0
            bit_llr = qam16_demod_bit_llr(log_post)
            code_llr = bit_llr[:n_code_bits][inv_perm]

            dec = conv_bcjr_decode(code_llr)
            hard = (dec.info_llr < 0).astype(np.uint8)
            bit_errors += int(np.sum(hard != info))
            total_bits += K
            if bit_errors >= min_errors:
                break
        ber = bit_errors / total_bits if total_bits > 0 else 0.0
        results.append({
            "eb_n0_db": eb_n0, "label": "no_isi",
            "bit_errors": bit_errors, "total_bits": total_bits, "ber": ber,
        })
        print(f"  [no-ISI] Eb/N0={eb_n0:.2f} dB  BER={ber:.6g}  "
              f"({bit_errors}/{total_bits})")
    return results


# ── Common result saving ────────────────────────────────────────

def _save_results(
    output_dir: Path, figure_name: str, eb_n0_range: list[float],
    labels: list[str], all_results: dict[str, dict[float, dict]],
    extra_config: dict,
) -> None:
    csv_path = output_dir / f"{figure_name}_results.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["eb_n0_db"] + labels)
        for snr in eb_n0_range:
            row = [snr]
            for label in labels:
                rec = all_results.get(label, {}).get(snr, {})
                row.append(f"{rec.get('ber', 0.0):.8g}")
            writer.writerow(row)

    np.savez(
        output_dir / f"{figure_name}_data.npz",
        eb_n0_db=np.array(eb_n0_range),
        **{label: np.array([all_results[label].get(s, {}).get("ber", 0.0)
                            for s in eb_n0_range])
           for label in labels},
    )

    with (output_dir / "config.json").open("w", encoding="utf-8") as f:
        json.dump(extra_config, f, indent=2)

    detail = {}
    for label in labels:
        detail[label] = [
            all_results[label].get(s)
            for s in eb_n0_range
            if s in all_results.get(label, {})
        ]
    with (output_dir / f"{figure_name}_detail.json").open("w", encoding="utf-8") as f:
        json.dump(detail, f, indent=2)

    print(f"\nResults saved to {output_dir}")


# ── Main ────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description="Parallel Monte Carlo runner")
    parser.add_argument(
        "figure", choices=["fig7", "fig8", "fig10", "fig11"],
    )
    parser.add_argument("--K", type=int, default=None)
    parser.add_argument("--frames", type=int, default=None)
    parser.add_argument("--procs", type=int, default=None)
    parser.add_argument("--seed", type=int, default=20260507)
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--m-states", type=int, default=None,
                        help="Override m_states for fig11")
    parser.add_argument("--future-len", type=int, default=None,
                        help="Override future_len for fig11")
    parser.add_argument("--isi-len", type=int, default=None,
                        help="Override isi_len for fig11")
    parser.add_argument("--turbo-iters", type=int, default=None,
                        help="Override turbo_iters for fig11")
    args = parser.parse_args()

    params = FIGURE_PARAMS[args.figure]
    eb_n0_range = list(params["eb_n0_range"])
    isi_len = params["isi_len"]
    turbo_iters = params["turbo_iters"]

    n_procs = args.procs or min(os.cpu_count() or 8, 52)
    output_dir = Path(args.output_dir or params["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    # ── Fig 7/8 ──
    if args.figure in ("fig7", "fig8"):
        tau = params["tau"]
        configs = params["configs"]
        K = args.K or 6000
        max_frames = args.frames or 100
        min_errors = 200 if K >= 3000 else max(50, K // 20)

        tasks = []
        task_id = 0
        for cfg in configs:
            for eb_n0 in eb_n0_range:
                tasks.append((
                    eb_n0, cfg["label"], cfg["m_states"], cfg["future_len"],
                    tau, isi_len, turbo_iters,
                    K, max_frames, min_errors,
                    args.seed + task_id,
                ))
                task_id += 1

        n_tasks = len(tasks)
        print(f"{'=' * 50}")
        print(f"Figure: {args.figure}  tau={tau}  isi_len={isi_len}")
        print(f"K={K}  frames<={max_frames}  min_errors={min_errors}")
        print(f"SNR range: Eb/N0 {eb_n0_range[0]:.1f}..{eb_n0_range[-1]:.1f} dB "
              f"(Es/N0 {eb_n0_range[0]+ES_N0_OFFSET_CONV:.1f}.."
              f"{eb_n0_range[-1]+ES_N0_OFFSET_CONV:.1f} dB)")
        print(f"Tasks: {n_tasks}  Workers: {n_procs}")
        print(f"{'=' * 50}\n")

        print("=== No-ISI baseline ===")
        t0 = time.time()
        noisi_results = _run_noisi_baseline_fig78(
            eb_n0_range, K, max_frames, min_errors, args.seed,
        )
        print(f"  Baseline done in {time.time() - t0:.1f}s\n")

        all_results: dict[str, dict[float, dict]] = {
            "no_isi": {r["eb_n0_db"]: r for r in noisi_results},
        }
        labels = ["no_isi"] + [c["label"] for c in configs]

        t_start = time.time()
        completed = 0
        with ProcessPoolExecutor(max_workers=n_procs) as pool:
            futures = {pool.submit(_run_one_task_fig78, task): task
                       for task in tasks}
            for future in as_completed(futures):
                result = future.result()
                label = result["label"]
                eb_n0 = result["eb_n0_db"]
                all_results.setdefault(label, {})[eb_n0] = result
                completed += 1
                elapsed = time.time() - t_start
                eta = elapsed / completed * (n_tasks - completed)
                print(f"  [{completed}/{n_tasks}] {label} Eb/N0={eb_n0:.2f} dB  "
                      f"BER={result['ber']:.6g}  "
                      f"({result['bit_errors']}/{result['total_bits']})  "
                      f"elapsed={elapsed:.0f}s  ETA={eta:.0f}s")

        print(f"\nAll FTN tasks done in {time.time() - t_start:.1f}s")

        _save_results(output_dir, args.figure, eb_n0_range, labels,
                      all_results, {
            "figure": args.figure, "tau": tau, "isi_len": isi_len,
            "turbo_iters": turbo_iters, "K": K, "frames": max_frames,
            "seed": args.seed, "n_procs": n_procs,
            "code_rate": CODE_RATE_CONV, "rolloff": 0.3,
            "pulse_span": 15, "sps": 128,
            "eb_n0_range": eb_n0_range, "configs": configs,
        })
        print("Done.")

    # ── Fig 10 ──
    elif args.figure == "fig10":
        tau_list = params["tau_list"]
        m_states = params["m_states"]
        future_len = params["future_len"]
        code_rate = params["code_rate"]

        K = args.K or 21842
        max_frames = args.frames or 50
        min_errors = 200 if K >= 5000 else max(50, K // 10)
        es_n0_offset = 10.0 * np.log10(code_rate)

        interleaver = s_random_interleaver(K, seed=args.seed + 100)

        configs = [{"label": f"tau={tau:.4f}", "tau": tau} for tau in tau_list]
        tasks = []
        task_id = 0
        for cfg in configs:
            for eb_n0 in eb_n0_range:
                tasks.append((
                    eb_n0, cfg["label"], cfg["tau"],
                    isi_len, m_states, future_len, turbo_iters,
                    K, max_frames, min_errors,
                    args.seed + task_id, code_rate,
                ))
                task_id += 1

        n_tasks = len(tasks)
        print(f"{'=' * 50}")
        print(f"Figure: fig10  Turbo R=1/3  isi_len={isi_len}  M={m_states}  L={future_len}")
        print(f"K={K}  frames<={max_frames}  turbo_iters={turbo_iters}  min_errors={min_errors}")
        print(f"SNR range: Eb/N0 {eb_n0_range[0]:.1f}..{eb_n0_range[-1]:.1f} dB "
              f"(Es/N0 {eb_n0_range[0]+es_n0_offset:.1f}.."
              f"{eb_n0_range[-1]+es_n0_offset:.1f} dB)")
        print(f"Tasks: {n_tasks}  Workers: {n_procs}")
        print(f"{'=' * 50}\n")

        print("=== No-ISI baseline ===")
        t0 = time.time()
        noisi_results = _run_noisi_baseline_fig10(
            eb_n0_range, K, max_frames, min_errors,
            args.seed, code_rate, interleaver,
        )
        print(f"  Baseline done in {time.time() - t0:.1f}s\n")

        all_results: dict[str, dict[float, dict]] = {
            "no_isi": {r["eb_n0_db"]: r for r in noisi_results},
        }
        labels = ["no_isi"] + [c["label"] for c in configs]

        t_start = time.time()
        completed = 0
        with ProcessPoolExecutor(max_workers=n_procs) as pool:
            futures = {pool.submit(_run_one_task_fig10, task): task
                       for task in tasks}
            for future in as_completed(futures):
                result = future.result()
                label = result["label"]
                eb_n0 = result["eb_n0_db"]
                all_results.setdefault(label, {})[eb_n0] = result
                completed += 1
                elapsed = time.time() - t_start
                eta = elapsed / completed * (n_tasks - completed)
                print(f"  [{completed}/{n_tasks}] {label} Eb/N0={eb_n0:.2f} dB  "
                      f"BER={result['ber']:.6g}  "
                      f"({result['bit_errors']}/{result['total_bits']})  "
                      f"elapsed={elapsed:.0f}s  ETA={eta:.0f}s")

        print(f"\nAll FTN tasks done in {time.time() - t_start:.1f}s")

        _save_results(output_dir, "fig10", eb_n0_range, labels,
                      all_results, {
            "figure": "fig10", "tau_list": tau_list, "isi_len": isi_len,
            "m_states": m_states, "future_len": future_len,
            "turbo_iters": turbo_iters, "code_rate": code_rate,
            "K": K, "frames": max_frames, "seed": args.seed,
            "n_procs": n_procs, "rolloff": 0.3, "pulse_span": 15, "sps": 128,
            "eb_n0_range": eb_n0_range,
            "interleaver": "random permutation of encoded bits before 16-QAM mapper",
            "interleaver_seed": "task seed + 1000",
        })
        print("Done.")

    # ── Fig 11 ──
    elif args.figure == "fig11":
        tau_list = params["tau_list"]
        m_states = args.m_states if args.m_states else params["m_states"]
        future_len = args.future_len if args.future_len else params["future_len"]
        isi_len = args.isi_len if args.isi_len else params["isi_len"]
        turbo_iters = args.turbo_iters if args.turbo_iters else params["turbo_iters"]
        code_rate = params["code_rate"]

        K = args.K or 6000
        max_frames = args.frames or 100
        min_errors = 200 if K >= 3000 else max(50, K // 10)
        es_n0_offset = 10.0 * np.log10(code_rate * 4)

        configs = [{"label": f"tau={tau:.4f}", "tau": tau} for tau in tau_list]
        tasks = []
        task_id = 0
        for cfg in configs:
            for eb_n0 in eb_n0_range:
                tasks.append((
                    eb_n0, cfg["label"], cfg["tau"],
                    isi_len, m_states, future_len, turbo_iters,
                    K, max_frames, min_errors,
                    args.seed + task_id, code_rate,
                ))
                task_id += 1

        n_tasks = len(tasks)
        print(f"{'=' * 50}")
        print(f"Figure: fig11  16-QAM + (7,5)  isi_len={isi_len}  M={m_states}  L={future_len}")
        print(f"K={K}  frames<={max_frames}  turbo_iters={turbo_iters}  min_errors={min_errors}")
        print(f"SNR range: Eb/N0 {eb_n0_range[0]:.1f}..{eb_n0_range[-1]:.1f} dB "
              f"(Es/N0 {eb_n0_range[0]+es_n0_offset:.1f}.."
              f"{eb_n0_range[-1]+es_n0_offset:.1f} dB)")
        print(f"Tasks: {n_tasks}  Workers: {n_procs}")
        print(f"{'=' * 50}\n")

        print("=== No-ISI baseline ===")
        t0 = time.time()
        noisi_results = _run_noisi_baseline_fig11(
            eb_n0_range, K, max_frames, min_errors,
            args.seed, code_rate,
        )
        print(f"  Baseline done in {time.time() - t0:.1f}s\n")

        all_results: dict[str, dict[float, dict]] = {
            "no_isi": {r["eb_n0_db"]: r for r in noisi_results},
        }
        labels = ["no_isi"] + [c["label"] for c in configs]

        t_start = time.time()
        completed = 0
        with ProcessPoolExecutor(max_workers=n_procs) as pool:
            futures = {pool.submit(_run_one_task_fig11, task): task
                       for task in tasks}
            for future in as_completed(futures):
                result = future.result()
                label = result["label"]
                eb_n0 = result["eb_n0_db"]
                all_results.setdefault(label, {})[eb_n0] = result
                completed += 1
                elapsed = time.time() - t_start
                eta = elapsed / completed * (n_tasks - completed)
                print(f"  [{completed}/{n_tasks}] {label} Eb/N0={eb_n0:.2f} dB  "
                      f"BER={result['ber']:.6g}  "
                      f"({result['bit_errors']}/{result['total_bits']})  "
                      f"elapsed={elapsed:.0f}s  ETA={eta:.0f}s")

        print(f"\nAll FTN tasks done in {time.time() - t_start:.1f}s")

        _save_results(output_dir, "fig11", eb_n0_range, labels,
                      all_results, {
            "figure": "fig11", "tau_list": tau_list, "isi_len": isi_len,
            "m_states": m_states, "future_len": future_len,
            "turbo_iters": turbo_iters, "code_rate": code_rate,
            "K": K, "frames": max_frames, "seed": args.seed,
            "n_procs": n_procs, "rolloff": 0.3, "pulse_span": 15, "sps": 128,
            "eb_n0_range": eb_n0_range,
        })
        print("Done.")


if __name__ == "__main__":
    main()
