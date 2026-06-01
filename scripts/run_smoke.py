from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys
from typing import Iterable

import numpy as np

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ftn.channel import ftn_awgn_channel
from ftn.equalizers.ungerboeck_mbcjr import ungerboeck_mbcjr_bpsk
from ftn.modulation import bpsk_hard_bits_from_llr, bpsk_modulate
from ftn.pulse import compute_g, generate_rrc


def _bpsk_from_bits(bits: np.ndarray) -> np.ndarray:
    return bpsk_modulate(bits)


def _hard_bits_from_llr(llr: np.ndarray) -> np.ndarray:
    return bpsk_hard_bits_from_llr(llr)


def run_uncoded_smoke(
    output_dir: str | Path = "results/smoke",
    seed: int = 20260506,
    frame_symbols: int = 200,
    frames_per_snr: int = 5,
    snr_db_values: Iterable[float] = (0.0, 2.0, 4.0, 6.0),
    tau: float = 0.8,
    rolloff: float = 0.3,
    pulse_span: int = 15,
    sps: int = 128,
    isi_len: int = 3,
    m_states: int = 8,
    future_len: int = 4,
) -> list[dict[str, float | int]]:
    """Run a deterministic short uncoded BPSK FTN smoke experiment."""
    if frame_symbols <= 0:
        raise ValueError("frame_symbols must be positive.")
    if frames_per_snr <= 0:
        raise ValueError("frames_per_snr must be positive.")

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(seed)
    t, h = generate_rrc(beta=rolloff, span=pulse_span, sps=sps)
    g = compute_g(t, h, tau=tau, isi_len=isi_len)

    rows: list[dict[str, float | int]] = []
    for snr_db in snr_db_values:
        n0 = 1.0 / (10.0 ** (float(snr_db) / 10.0))
        bit_errors = 0
        total_bits = 0
        for _ in range(frames_per_snr):
            bits = rng.integers(0, 2, size=frame_symbols, dtype=np.uint8)
            symbols = _bpsk_from_bits(bits)
            y = ftn_awgn_channel(symbols, g, n0=n0, rng=rng)
            result = ungerboeck_mbcjr_bpsk(
                y,
                g,
                n0=n0,
                isi_len=isi_len,
                m_states=m_states,
                future_len=future_len,
                initial_state=tuple(1.0 for _ in range(isi_len)),
            )
            hard = _hard_bits_from_llr(result.llr)
            bit_errors += int(np.count_nonzero(hard != bits))
            total_bits += int(bits.size)

        rows.append(
            {
                "snr_db": float(snr_db),
                "bit_errors": bit_errors,
                "total_bits": total_bits,
                "ber": bit_errors / float(total_bits),
            }
        )

    csv_path = output_path / "uncoded_smoke.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["snr_db", "bit_errors", "total_bits", "ber"])
        writer.writeheader()
        writer.writerows(rows)

    config = {
        "seed": seed,
        "frame_symbols": frame_symbols,
        "frames_per_snr": frames_per_snr,
        "snr_db_values": [float(v) for v in snr_db_values],
        "tau": tau,
        "rolloff": rolloff,
        "pulse_span": pulse_span,
        "sps": sps,
        "isi_len": isi_len,
        "m_states": m_states,
        "future_len": future_len,
        "g": {str(k): float(v) for k, v in sorted(g.items())},
    }
    with (output_path / "uncoded_smoke_config.json").open("w", encoding="utf-8") as f:
        json.dump(config, f, indent=2)

    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Run a short uncoded FTN M-BCJR smoke test.")
    parser.add_argument("--output-dir", default="results/smoke")
    parser.add_argument("--seed", type=int, default=20260506)
    parser.add_argument("--frame-symbols", type=int, default=200)
    parser.add_argument("--frames-per-snr", type=int, default=5)
    parser.add_argument("--tau", type=float, default=0.8)
    parser.add_argument("--isi-len", type=int, default=3)
    parser.add_argument("--m-states", type=int, default=8)
    parser.add_argument("--future-len", type=int, default=4)
    parser.add_argument(
        "--snr-db",
        type=float,
        nargs="+",
        default=[0.0, 2.0, 4.0, 6.0],
        help="SNR points in dB.",
    )
    args = parser.parse_args()

    rows = run_uncoded_smoke(
        output_dir=args.output_dir,
        seed=args.seed,
        frame_symbols=args.frame_symbols,
        frames_per_snr=args.frames_per_snr,
        snr_db_values=args.snr_db,
        tau=args.tau,
        isi_len=args.isi_len,
        m_states=args.m_states,
        future_len=args.future_len,
    )
    for row in rows:
        print(
            f"SNR={row['snr_db']:.2f} dB "
            f"errors={row['bit_errors']}/{row['total_bits']} "
            f"BER={row['ber']:.6g}"
        )


if __name__ == "__main__":
    main()
