"""Diagnostic: test fig11 τ=2/3 with different parameter combinations.

Tests a single SNR point (12 dB) with K=3000, 20 frames, varying:
- m_states: 4, 8, 16
- future_len: 3, 4
- llr_clip: 8, 20
- turbo_iters: 10, 30
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ftn.channel import ftn_awgn_channel_complex
from ftn.coding.conv import conv_bcjr_decode, conv_encode_75
from ftn.equalizers.simplified_mbcjr import simplified_mbcjr_qam16
from ftn.modulation import qam16_modulate
from ftn.pulse import compute_g, generate_rrc


def run_one_config(
    tau: float, eb_n0_db: float, K: int, n_frames: int,
    m_states: int, future_len: int, turbo_iters: int, llr_clip: float,
    isi_len: int, seed: int,
) -> dict:
    code_rate = 0.5
    es_n0_offset = 10.0 * np.log10(code_rate * 4)
    es_n0 = eb_n0_db + es_n0_offset
    n0 = 1.0 / (10.0 ** (es_n0 / 10.0))

    t, h = generate_rrc(beta=0.3, span=15, sps=128)
    g = compute_g(t, h, tau=tau, isi_len=isi_len)
    rng = np.random.default_rng(seed)

    bit_errors = 0
    total_bits = 0
    t0 = time.time()

    for fi in range(n_frames):
        info_bits = rng.integers(0, 2, size=K, dtype=np.uint8)
        encoded = conv_encode_75(info_bits)
        n_code_bits = encoded.size
        pad = (4 - n_code_bits % 4) % 4
        if pad > 0:
            tx_bits = np.concatenate([encoded, np.zeros(pad, dtype=np.uint8)])
        else:
            tx_bits = encoded
        n_sym = tx_bits.size // 4

        symbols = qam16_modulate(tx_bits)
        y = ftn_awgn_channel_complex(symbols, g, n0=n0, rng=rng)

        la = np.zeros(n_sym * 4, dtype=np.float64)

        for _it in range(turbo_iters + 1):
            det_result = simplified_mbcjr_qam16(
                y, g, n0=n0, la=la,
                isi_len=isi_len, m_states=m_states, future_len=future_len,
            )
            det_llr = np.clip(det_result.bit_llr, -llr_clip, llr_clip)
            det_ext = np.clip(det_llr - la, -llr_clip, llr_clip)
            code_llr_in = det_ext[:n_code_bits]

            dec = conv_bcjr_decode(code_llr_in)
            info_llr = np.clip(dec.info_llr, -llr_clip, llr_clip)
            code_llr = np.clip(dec.code_llr, -llr_clip, llr_clip)
            dec_ext = np.clip(code_llr - code_llr_in, -llr_clip, llr_clip)

            la = np.zeros(n_sym * 4, dtype=np.float64)
            la[:n_code_bits] = dec_ext

        hard_info = (info_llr < 0.0).astype(np.uint8)
        bit_errors += int(np.sum(hard_info != info_bits))
        total_bits += int(info_bits.size)

    elapsed = time.time() - t0
    ber = bit_errors / total_bits if total_bits > 0 else 0.0
    return {
        "ber": ber, "errors": bit_errors, "bits": total_bits,
        "elapsed": elapsed,
    }


def main():
    K = 3000
    n_frames = 20
    eb_n0_db = 12.0
    tau = 2 / 3
    isi_len = 3
    seed = 20260523

    configs = [
        # Label, m_states, future_len, turbo_iters, llr_clip
        ("M4_L3_clip20_it10",   4,  3, 10, 20.0),
        ("M4_L3_clip8_it10",    4,  3, 10,  8.0),
        ("M4_L3_clip20_it30",   4,  3, 30, 20.0),
        ("M4_L3_clip8_it30",    4,  3, 30,  8.0),
        ("M8_L3_clip20_it10",   8,  3, 10, 20.0),
        ("M8_L3_clip8_it10",    8,  3, 10,  8.0),
        ("M16_L3_clip20_it10", 16,  3, 10, 20.0),
        ("M16_L3_clip8_it10",  16,  3, 10,  8.0),
        ("M4_L4_clip20_it10",   4,  4, 10, 20.0),
        ("M4_L4_clip8_it10",    4,  4, 10,  8.0),
        ("M8_L4_clip8_it10",    8,  4, 10,  8.0),
    ]

    print(f"Fig 11 diagnostic: tau=2/3, Eb/N0={eb_n0_db} dB, K={K}, frames={n_frames}")
    print(f"{'Config':<25} {'BER':>10} {'Errors':>8} {'Bits':>8} {'Time':>8}")
    print("-" * 70)

    for label, m_states, future_len, turbo_iters, llr_clip in configs:
        result = run_one_config(
            tau, eb_n0_db, K, n_frames,
            m_states, future_len, turbo_iters, llr_clip,
            isi_len, seed,
        )
        print(f"{label:<25} {result['ber']:10.6f} {result['errors']:8d} "
              f"{result['bits']:8d} {result['elapsed']:7.1f}s")

    print("\nDone.")


if __name__ == "__main__":
    main()
