"""Generate BER vs Eb/N0 plots from merged CSV results.

Usage:
    python scripts/plot_results.py fig7
    python scripts/plot_results.py fig8
    python scripts/plot_results.py both
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

PROJECT = Path(__file__).resolve().parent.parent

MARKERS = {
    "no_isi": {"marker": "o", "color": "black", "linestyle": "-", "label": "No-ISI (conv. code)"},
    "M2_L3": {"marker": "s", "color": "blue", "linestyle": "--", "label": "M-BCJR M=2, L=3"},
    "M2_L5": {"marker": "^", "color": "red", "linestyle": "-.", "label": "M-BCJR M=2, L=5"},
    "M8_L5": {"marker": "D", "color": "blue", "linestyle": "--", "label": "M-BCJR M=8, L=5"},
    "M8_L7": {"marker": "v", "color": "red", "linestyle": "-.", "label": "M-BCJR M=8, L=7"},
    "no_isi": {"marker": "o", "color": "black", "linestyle": "-", "label": "No-ISI (16-QAM)"},
    "tau06667": {"marker": "s", "color": "blue", "linestyle": "--", "label": "FTN τ=2/3, M=4, L=3"},
    "tau08000": {"marker": "^", "color": "red", "linestyle": "-.", "label": "FTN τ=0.8, M=4, L=3"},
    "tau10000": {"marker": "s", "color": "blue", "linestyle": "--", "label": "FTN τ=1.0, M=8, L=5"},
    "tau06667_turbo": {"marker": "D", "color": "green", "linestyle": "-.", "label": "FTN τ=2/3, M=8, L=5"},
    "tau05000": {"marker": "v", "color": "red", "linestyle": ":", "label": "FTN τ=0.5, M=8, L=5"},
}

FIG_CONFIG = {
    "fig7": {
        "result_dir": "results/fig7_tau_05",
        "csv_name": "fig7_results.csv",
        "png_name": "fig7_ber.png",
        "title": "Fig. 7 — BPSK, τ=0.5, (7,5) conv. code, turbo iters=5",
        "curves": ["no_isi", "M2_L3", "M2_L5"],
    },
    "fig8": {
        "result_dir": "results/fig8_tau_035",
        "csv_name": "fig8_results.csv",
        "png_name": "fig8_ber.png",
        "title": "Fig. 8 — BPSK, τ=0.35, (7,5) conv. code, turbo iters=15",
        "curves": ["no_isi", "M8_L5", "M8_L7"],
    },
    "fig10": {
        "result_dir": "results/fig10",
        "csv_name": "fig10_results.csv",
        "png_name": "fig10_ber.png",
        "title": "Fig. 10 — BPSK, turbo R=1/3, turbo equalization",
        "curves": ["no_isi", "tau10000", "tau06667", "tau05000"],
    },
    "fig11": {
        "result_dir": "results/fig11",
        "csv_name": "fig11_results.csv",
        "png_name": "fig11_ber.png",
        "title": "Fig. 11 — 16-QAM, simplified M-BCJR, M=4, L=3",
        "curves": ["no_isi", "tau06667", "tau08000"],
    },
}


def load_csv(path: Path) -> dict[str, np.ndarray]:
    data = np.genfromtxt(path, delimiter=",", names=True)
    result = {"eb_n0_db": data["eb_n0_db"]}
    for name in data.dtype.names[1:]:
        result[name] = data[name]
    return result


def plot_figure(fig_key: str) -> None:
    cfg = FIG_CONFIG[fig_key]
    result_dir = PROJECT / cfg["result_dir"]
    csv_path = result_dir / cfg["csv_name"]
    png_path = result_dir / cfg["png_name"]

    if not csv_path.exists():
        print(f"CSV not found: {csv_path}")
        return

    data = load_csv(csv_path)
    snr = data["eb_n0_db"]

    fig, ax = plt.subplots(figsize=(8, 5.5))

    for curve_key in cfg["curves"]:
        if curve_key not in data:
            print(f"  Warning: curve '{curve_key}' not in CSV, skipping")
            continue
        ber = data[curve_key]
        style = MARKERS[curve_key]
        ax.semilogy(snr, ber, marker=style["marker"], color=style["color"],
                    linestyle=style["linestyle"], linewidth=1.5, markersize=6,
                    label=style["label"])

    ax.set_xlabel("$E_b/N_0$ (dB)", fontsize=12)
    ax.set_ylabel("BER", fontsize=12)
    ax.set_title(cfg["title"], fontsize=12)
    ax.set_xlim(snr[0] - 0.05, snr[-1] + 0.05)
    ax.set_ylim(1e-6, 1)
    ax.grid(True, which="both", linestyle=":", alpha=0.5)
    ax.legend(fontsize=10, loc="upper right")
    fig.tight_layout()

    fig.savefig(png_path, dpi=200)
    print(f"Saved: {png_path}")
    plt.close(fig)


def main() -> None:
    which = sys.argv[1] if len(sys.argv) > 1 else "both"
    if which == "both":
        for k in FIG_CONFIG:
            plot_figure(k)
    elif which in FIG_CONFIG:
        plot_figure(which)
    else:
        print(f"Unknown figure '{which}'. Use: fig7, fig8, or both")
        sys.exit(1)


if __name__ == "__main__":
    main()
