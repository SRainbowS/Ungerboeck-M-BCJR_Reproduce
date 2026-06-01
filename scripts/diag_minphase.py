"""Diagnostic: fix min_phase_from_pulse for tau=0.35."""
import sys
import numpy as np
from pathlib import Path

PROJECT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT))

from ftn.baselines.forney_model import (
    min_phase_from_pulse, sample_pulse_at_ftn_rate, spectral_factorize
)
from ftn.pulse import compute_g, generate_rrc

for tau in [0.5, 0.35]:
    isi_len = 7 if tau == 0.5 else 10
    t, h = generate_rrc(beta=0.3, span=15, sps=128)
    g = compute_g(t, h, tau=tau, isi_len=isi_len)
    print(f"\n=== tau={tau}, isi_len={isi_len} ===")
    print(f"g[0] = {g[0]:.6f}")

    # Method 1: current (hilbert cepstral)
    c = sample_pulse_at_ftn_rate(t, h, tau=tau)
    print(f"Sampled pulse c: length={len(c)}, c[center]={c[len(c)//2]:.6f}")
    print(f"c has_nan={np.any(np.isnan(c))}, sum_c={np.sum(c):.6f}")

    # Check spectrum of c
    C_w = np.fft.rfft(c)
    print(f"C_w: min_abs={np.min(np.abs(C_w)):.6e}, has_zero={np.any(np.abs(C_w) < 1e-10)}")

    # Try different methods
    from scipy.signal import minimum_phase
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        try:
            v_hilbert = minimum_phase(c, method="hilbert")
            v_hilbert = v_hilbert * np.sqrt(g[0] / np.sum(v_hilbert**2)) if np.sum(v_hilbert**2) > 0 else v_hilbert
            print(f"  hilbert: v[0]={v_hilbert[0]:.6f}, len={len(v_hilbert)}, has_nan={np.any(np.isnan(v_hilbert))}")
        except Exception as e:
            print(f"  hilbert: FAILED - {e}")

        try:
            v_homomorphic = minimum_phase(c, method="homomorphic")
            v_homomorphic = v_homomorphic * np.sqrt(g[0] / np.sum(v_homomorphic**2)) if np.sum(v_homomorphic**2) > 0 else v_homomorphic
            print(f"  homomorphic: v[0]={v_homomorphic[0]:.6f}, len={len(v_homomorphic)}, has_nan={np.any(np.isnan(v_homomorphic))}")
        except Exception as e:
            print(f"  homomorphic: FAILED - {e}")

    # Method 2: spectral factorization from g-taps
    v_spec = spectral_factorize(g)
    print(f"  spectral_factorize: v[0]={v_spec[0]:.6f}, len={len(v_spec)}, has_nan={np.any(np.isnan(v_spec))}")
    print(f"  sum(v^2) = {np.sum(v_spec**2):.6f}")

    # Method 3: root-based minimum phase from g-taps
    # Build polynomial from g, find roots, select those inside unit circle
    L = isi_len
    # G(z) = sum_{k=-L}^{L} g[k] z^{-k} = z^{-L} * P(z) where P(z) = sum g[k] z^{L-k}
    poly = np.zeros(2*L+1)
    for k, val in g.items():
        poly[L - k] = val
    roots = np.roots(poly)
    # Select L roots inside or on unit circle
    inside = roots[np.abs(roots) <= 1.0 + 1e-6]
    if len(inside) > L:
        # Pick L closest to origin
        inside = inside[np.argsort(np.abs(inside))[:L]]
    if len(inside) < L:
        print(f"  root-based: only {len(inside)}/{L} roots inside, skipping")
    else:
        v_poly = np.real(np.poly(inside))
        # Scale: g[0] = sum(v^2)
        scale = np.sqrt(g[0] / np.sum(v_poly**2))
        v_root = v_poly * scale
        if v_root[0] < 0:
            v_root = -v_root
        print(f"  root-based: v[0]={v_root[0]:.6f}, len={len(v_root)}, has_nan={np.any(np.isnan(v_root))}")
        print(f"  v_root[:5] = {v_root[:5]}")

print("\nDone.")
