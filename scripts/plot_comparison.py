"""Plot overlay BER curves for algorithm comparison.

Usage:
    python scripts/plot_comparison.py uncoded --tau 0.5
    python scripts/plot_comparison.py uncoded --tau 0.35
    python scripts/plot_comparison.py turbo --tau 0.5
    python scripts/plot_comparison.py turbo --tau 0.35
    python scripts/plot_comparison.py all
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

PROJECT = Path(__file__).resolve().parent.parent

# ---------------------------------------------------------------------------
# Style definitions for each curve
# ---------------------------------------------------------------------------

UNCODED_STYLES = {
    "full_bcjr":  {"marker": "o",  "color": "black",   "ls": "-",  "label": "Full BCJR (oracle)"},
    "ung_M4":     {"marker": "s",  "color": "#1f77b4",  "ls": "--", "label": "Ungerboeck M-BCJR M=4"},
    "ung_M8":     {"marker": "D",  "color": "#1f77b4",  "ls": "-",  "label": "Ungerboeck M-BCJR M=8"},
    "prlja_M4":   {"marker": "^",  "color": "#d62728",  "ls": "--", "label": "Paper [14] M-BCJR M=4"},
    "prlja_M8":   {"marker": "v",  "color": "#d62728",  "ls": "-",  "label": "Paper [14] M-BCJR M=8"},
    "cs_nu2":     {"marker": "P",  "color": "#2ca02c",  "ls": "--", "label": "Paper [26] Shortened nu=2"},
    "cs_nu3":     {"marker": "X",  "color": "#2ca02c",  "ls": "-.", "label": "Paper [26] Shortened nu=3"},
    "cs_nu5":     {"marker": "*",  "color": "#2ca02c",  "ls": "-",  "label": "Paper [26] Shortened nu=5"},
}

TURBO_STYLES = {
    "no_isi":     {"marker": "o",  "color": "black",   "ls": "-",  "label": "No-ISI baseline"},
    "ung_M4":     {"marker": "s",  "color": "#1f77b4",  "ls": "--", "label": "Ungerboeck M-BCJR M=4"},
    "ung_M8":     {"marker": "D",  "color": "#1f77b4",  "ls": "-",  "label": "Ungerboeck M-BCJR M=8"},
    "prlja_M4":   {"marker": "^",  "color": "#d62728",  "ls": "--", "label": "Paper [14] Backup M-BCJR M=4"},
    "prlja_M8":   {"marker": "v",  "color": "#d62728",  "ls": "-",  "label": "Paper [14] Backup M-BCJR M=8"},
    "cs_nu2":     {"marker": "P",  "color": "#2ca02c",  "ls": "--", "label": "Paper [26] Shortened nu=2"},
    "cs_nu3":     {"marker": "X",  "color": "#2ca02c",  "ls": "-.", "label": "Paper [26] Shortened nu=3"},
}


def _load_csv(path: Path) -> dict[str, np.ndarray]:
    data = np.genfromtxt(path, delimiter=",", names=True)
    result = {}
    for name in data.dtype.names:
        result[name] = data[name]
    return result


def _plot_one(data: dict, styles: dict, title: str, png_path: Path,
              snr_key: str = "snr_db"):
    fig, ax = plt.subplots(figsize=(9, 6))

    snr = data[snr_key]
    for key, style in styles.items():
        if key not in data:
            print(f"  Warning: '{key}' not in data, skipping")
            continue
        ber = data[key]
        mask = ber > 0
        if mask.sum() == 0:
            continue
        ax.semilogy(snr[mask], ber[mask],
                    marker=style["marker"], color=style["color"],
                    linestyle=style["ls"], linewidth=1.5, markersize=6,
                    label=style["label"])

    ax.set_xlabel("$E_b/N_0$ (dB)", fontsize=12)
    ax.set_ylabel("BER", fontsize=12)
    ax.set_title(title, fontsize=12)
    ax.set_xlim(snr[0] - 0.05, snr[-1] + 0.05)
    ax.set_ylim(1e-6, 1)
    ax.grid(True, which="both", linestyle=":", alpha=0.5)
    ax.legend(fontsize=9, loc="upper right")
    fig.tight_layout()
    fig.savefig(png_path, dpi=200)
    print(f"Saved: {png_path}")
    plt.close(fig)


def main():
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    tau = 0.5
    for i, arg in enumerate(sys.argv[1:], 1):
        if arg == "--tau" and i < len(sys.argv) - 1:
            tau = float(sys.argv[i + 1])

    tau_tag = f"tau{int(tau * 1000):04d}"

    if what in ("uncoded", "all"):
        csv_path = PROJECT / f"results/comparison/uncoded_{tau_tag}" / "uncoded_ber.csv"
        if csv_path.exists():
            data = _load_csv(csv_path)
            _plot_one(data, UNCODED_STYLES,
                      f"Uncoded BPSK BER Comparison — FTN tau={tau}",
                      csv_path.parent / "uncoded_ber.png")
        else:
            print(f"Not found: {csv_path}")

    if what in ("turbo", "all"):
        csv_path = PROJECT / f"results/comparison/turbo_{tau_tag}" / "turbo_ber.csv"
        if csv_path.exists():
            data = _load_csv(csv_path)
            _plot_one(data, TURBO_STYLES,
                      f"Turbo Equalization BER — FTN tau={tau}, (7,5) conv",
                      csv_path.parent / "turbo_ber.png")
        else:
            print(f"Not found: {csv_path}")


if __name__ == "__main__":
    main()
