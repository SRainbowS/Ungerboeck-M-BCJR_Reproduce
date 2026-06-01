"""Merge partial Monte Carlo results from parallel workers into final CSV/NPZ.

Usage:
    python scripts/merge_results.py <partial_dir> <output_dir> <fig7|fig8>
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

# Ensure project root is on sys.path so `ftn` can be imported
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np


def main() -> None:
    if len(sys.argv) < 4:
        print("Usage: merge_results.py <partial_dir> <output_dir> <fig7|fig8>")
        sys.exit(1)

    partial_dir = Path(sys.argv[1])
    output_dir = Path(sys.argv[2])
    figure = sys.argv[3]

    if not partial_dir.exists():
        print(f"Partial directory not found: {partial_dir}")
        sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)

    # Read all partial CSV files
    records: dict[str, dict[float, dict]] = {}
    for csv_file in sorted(partial_dir.glob("*.csv")):
        with csv_file.open("r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                label = row["label"]
                snr = float(row["eb_n0_db"])
                if label not in records:
                    records[label] = {}
                records[label][snr] = {
                    "eb_n0_db": snr,
                    "label": label,
                    "m_states": int(row["m_states"]),
                    "future_len": int(row["future_len"]),
                    "bit_errors": int(row["bit_errors"]),
                    "total_bits": int(row["total_bits"]),
                    "ber": float(row["ber"]),
                }

    if not records:
        print("No partial results found.")
        sys.exit(1)

    # Collect all SNR points across all curves
    all_snr = sorted(set(
        snr for label_data in records.values() for snr in label_data
    ))

    # Also generate no-ISI baseline (fast, do it here)
    print("Computing no-ISI baseline...")
    from ftn.coding.conv import conv_bcjr_decode, conv_encode_75
    from ftn.modulation import bpsk_modulate

    K = 6000
    frames_per_snr = 100
    rng = np.random.default_rng(20260507)
    code_rate = 0.5

    noisi_records = {}
    for snr in all_snr:
        es_n0 = snr + 10 * np.log10(code_rate)
        n0 = 1.0 / (10 ** (es_n0 / 10.0))
        bit_errors = 0
        total_bits = 0
        for _ in range(frames_per_snr):
            info = rng.integers(0, 2, size=K, dtype=np.uint8)
            enc = conv_encode_75(info)
            sym = bpsk_modulate(enc)
            y = sym + rng.standard_normal(sym.size) * np.sqrt(n0 / 2.0)
            code_llr = 4.0 * y / n0
            dec = conv_bcjr_decode(code_llr)
            hard = (dec.info_llr < 0).astype(np.uint8)
            bit_errors += int(np.sum(hard != info))
            total_bits += K
            if bit_errors >= 200 and total_bits >= 5 * K:
                break
        ber = bit_errors / total_bits if total_bits > 0 else 0.0
        noisi_records[snr] = {
            "eb_n0_db": snr, "label": "no_isi",
            "bit_errors": bit_errors, "total_bits": total_bits, "ber": ber,
        }
        print(f"  [no-ISI] Eb/N0={snr:.2f} dB  BER={ber:.6g}  ({bit_errors}/{total_bits})")

    records["no_isi"] = noisi_records

    # Build ordered label list
    labels = ["no_isi"] + sorted(l for l in records if l != "no_isi")

    # ---- Save CSV ----
    csv_path = output_dir / f"{figure}_results.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["eb_n0_db"] + labels)
        for snr in all_snr:
            row = [snr]
            for label in labels:
                rec = records.get(label, {}).get(snr, {})
                row.append(f"{rec.get('ber', 0.0):.8g}")
            writer.writerow(row)
    print(f"\nSaved CSV to {csv_path}")

    # ---- Save NPZ ----
    np.savez(
        output_dir / f"{figure}_data.npz",
        eb_n0_db=np.array(all_snr),
        **{label: np.array([records[label].get(snr, {}).get("ber", 0.0) for snr in all_snr])
           for label in labels},
    )

    # ---- Save detailed JSON ----
    detail = {}
    for label in labels:
        detail[label] = [
            records[label][snr] for snr in all_snr if snr in records.get(label, {})
        ]
    with (output_dir / f"{figure}_detail.json").open("w", encoding="utf-8") as f:
        json.dump(detail, f, indent=2)

    # ---- Save config ----
    fig_params = {
        "fig7": {"tau": 0.5, "isi_len": 7, "turbo_iters": 5, "configs": [
            {"label": "M2_L3", "m_states": 2, "future_len": 3},
            {"label": "M2_L5", "m_states": 2, "future_len": 5},
        ]},
        "fig8": {"tau": 0.35, "isi_len": 10, "turbo_iters": 15, "configs": [
            {"label": "M8_L5", "m_states": 8, "future_len": 5},
            {"label": "M8_L7", "m_states": 8, "future_len": 7},
        ]},
    }[figure]

    config = {
        "figure": figure,
        "K": K, "code_rate": 0.5, "rolloff": 0.3, "pulse_span": 15, "sps": 128,
        "llr_clip": 20.0, "snr_range": all_snr,
        "eb_n0_offset_db": float(10 * np.log10(code_rate)),
        **fig_params,
    }
    with (output_dir / "config.json").open("w", encoding="utf-8") as f:
        json.dump(config, f, indent=2)

    print(f"Results saved to {output_dir}")


if __name__ == "__main__":
    main()
