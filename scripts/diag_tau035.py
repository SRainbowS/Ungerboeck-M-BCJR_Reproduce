"""Diagnostic: why does CS equalizer give BER=1.0 for tau=0.35?"""
import sys
import numpy as np
from pathlib import Path

PROJECT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT))

from ftn.baselines.forney_model import forney_channel, min_phase_from_pulse
from ftn.baselines.channel_shortening import compute_shortened_params_from_v, apply_shortening_filter_fft
from ftn.baselines.shortened_bcjr import cs_equalizer_bpsk
from ftn.channel import ftn_awgn_channel
from ftn.equalizers.full_bcjr import full_bcjr_bpsk
from ftn.pulse import compute_g, generate_rrc

for tau in [0.5, 0.35]:
    isi_len = 7 if tau == 0.5 else 10
    t, h = generate_rrc(beta=0.3, span=15, sps=128)
    g = compute_g(t, h, tau=tau, isi_len=isi_len)
    v = min_phase_from_pulse(t, h, tau=tau, g0=g[0])

    print(f"\n=== tau={tau}, isi_len={isi_len} ===")
    print(f"v length: {len(v)}")
    print(f"v[0] = {v[0]:.6f}")
    print(f"v[:8] = {v[:8]}")
    print(f"sum(v^2) = {np.sum(v**2):.6f}, g[0] = {g[0]:.6f}")

    n_bits = 500
    rng = np.random.default_rng(42)
    x = rng.choice([-1.0, 1.0], size=n_bits)
    snr_db = 5.0
    n0 = 1.0 / (10.0 ** (snr_db / 10.0))

    # Forney channel
    y = forney_channel(x, v, n0=n0, rng=rng)
    print(f"y stats: mean={y.mean():.4f}, std={y.std():.4f}")

    # Test CS for different nu
    for nu in [2, 3, 5]:
        if nu >= isi_len:
            continue
        n_fft = 1
        while n_fft < n_bits + len(v):
            n_fft *= 2
        Z_w, g_diag, g_off, n_fft = compute_shortened_params_from_v(v, n0, nu, n_fft)

        z = apply_shortening_filter_fft(y, Z_w, n_fft)
        print(f"\n  CS nu={nu}:")
        print(f"    g_diag = {g_diag:.6f}")
        print(f"    g_off = {g_off}")
        print(f"    Z_w: has_nan={np.any(np.isnan(Z_w))}, has_inf={np.any(np.isinf(Z_w))}, "
              f"min={np.nanmin(Z_w):.4f}, max={np.nanmax(Z_w):.4f}")
        print(f"    z: has_nan={np.any(np.isnan(z))}, has_inf={np.any(np.isinf(z))}, "
              f"mean={np.nanmean(z):.4f}, std={np.nanstd(z):.4f}")

        # Run CS equalizer
        result = cs_equalizer_bpsk(y, v, n0, nu)
        decisions = np.sign(result.llr)
        ber = np.mean(decisions != x)
        print(f"    BER = {ber:.4f}")

        # Check if flipping LLR fixes it
        ber_flip = np.mean(-decisions != x)
        print(f"    BER (flipped) = {ber_flip:.4f}")

    # Full BCJR for reference
    result_full = full_bcjr_bpsk(y, g, n0=n0, isi_len=isi_len)
    ber_full = np.mean(np.sign(result_full.llr) != x)
    print(f"\n  Full BCJR BER = {ber_full:.4f}")

print("\nDone.")
