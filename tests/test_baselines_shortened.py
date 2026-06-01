"""Tests for Paper [26] channel shortening baseline algorithms."""

import pytest
import numpy as np

from ftn.baselines.channel_shortening import (
    compute_b_matrix,
    solve_banded_g,
    compute_shortened_params,
    apply_shortening_filter,
    compute_shortened_params_from_v,
    apply_shortening_filter_fft,
)
from ftn.baselines.shortened_bcjr import (
    shortened_bcjr_bpsk,
    cs_equalizer_bpsk,
)
from ftn.channel import build_isi_matrix


# ---------------------------------------------------------------------------
# B matrix tests
# ---------------------------------------------------------------------------

def test_b_matrix_properties():
    """B should be symmetric and positive semi-definite."""
    g = {0: 1.0, 1: 0.3, -1: 0.3, 2: 0.1, -2: 0.1}
    B = compute_b_matrix(g, n0=1.0, block_len=50)
    np.testing.assert_allclose(B, B.T, atol=1e-10)
    eigvals = np.linalg.eigvalsh(B)
    assert np.all(eigvals >= -1e-8)


def test_b_matrix_noise_limits():
    """n0 -> 0: B -> 0 (MMSE vanishes); n0 -> inf: B -> I."""
    g = {0: 1.0, 1: 0.3, -1: 0.3}
    B_low = compute_b_matrix(g, n0=0.001, block_len=30)
    assert np.linalg.norm(B_low) < 1.0
    B_high = compute_b_matrix(g, n0=1000.0, block_len=30)
    np.testing.assert_allclose(B_high, np.eye(30), atol=0.1)


# ---------------------------------------------------------------------------
# G matrix tests
# ---------------------------------------------------------------------------

def test_g_banded_property():
    """G should have bandwidth nu (entries outside band are zero)."""
    g = {0: 1.0, 1: 0.3, -1: 0.3, 2: 0.1, -2: 0.1}
    B = compute_b_matrix(g, n0=0.5, block_len=30)
    nu = 2
    G = solve_banded_g(B, nu)
    N = G.shape[0]
    for i in range(N):
        for j in range(N):
            if abs(i - j) > nu:
                assert abs(G[i, j]) < 1e-10


def test_g_positive_definite():
    """G should be positive definite."""
    g = {0: 1.0, 1: 0.3, -1: 0.3, 2: 0.1, -2: 0.1}
    B = compute_b_matrix(g, n0=0.5, block_len=30)
    G = solve_banded_g(B, nu=2)
    eigvals = np.linalg.eigvalsh(G)
    assert np.all(eigvals > 0)


def test_gmi_monotonic_with_nu():
    """GMI should be non-decreasing with bandwidth nu."""
    g = {0: 1.0, 1: 0.3, -1: 0.3, 2: 0.1, -2: 0.1}
    B = compute_b_matrix(g, n0=0.5, block_len=30)
    gmis = []
    for nu in [1, 2, 3, 4]:
        G = solve_banded_g(B, nu)
        sign, logdet = np.linalg.slogdet(G)
        gmi = logdet - np.trace(G @ B)
        gmis.append(gmi)
    for i in range(1, len(gmis)):
        assert gmis[i] >= gmis[i - 1] - 1e-4


# ---------------------------------------------------------------------------
# Shortened BCJR tests
# ---------------------------------------------------------------------------

def test_shortened_bcjr_finite_llr():
    """LLRs from shortened BCJR should all be finite."""
    g = {0: 1.0, 1: 0.3, -1: 0.3, 2: 0.1, -2: 0.1}
    h_r, g_diag, g_off = compute_shortened_params(g, n0=0.5, nu=2)
    rng = np.random.default_rng(42)
    y = rng.standard_normal(50)
    z = apply_shortening_filter(y, h_r)
    result = shortened_bcjr_bpsk(z, g_diag, g_off, nu=2)
    assert np.all(np.isfinite(result.llr))


def test_shortened_bcjr_ber_with_snr():
    """BER should decrease with decreasing noise (increasing SNR)."""
    g = {0: 1.0, 1: 0.3, -1: 0.3}
    bers = []
    for n0 in [2.0, 0.5, 0.1]:
        errors = 0
        total = 0
        rng = np.random.default_rng(42)
        h_r, g_diag, g_off = compute_shortened_params(g, n0, nu=2)
        for _ in range(10):
            x = rng.choice([-1, 1], size=200)
            H = build_isi_matrix(g, len(x))
            n_vec = np.sqrt(n0) * rng.standard_normal(len(x))
            y = H @ x + n_vec
            z = apply_shortening_filter(y, h_r)
            result = shortened_bcjr_bpsk(z, g_diag, g_off, nu=2)
            errors += np.sum(np.sign(result.llr) != x)
            total += len(x)
        bers.append(errors / total)
    assert bers[-1] < bers[0]


def test_cs_equalizer_pipeline():
    """Full cs_equalizer_bpsk pipeline (Forney model) should produce finite
    LLRs and detect more than half the symbols correctly at moderate SNR."""
    # Use a simple 3-tap min-phase channel (causal)
    v = np.array([1.0, 0.3, 0.1])
    v = v / np.sqrt(np.sum(v**2))  # normalise energy to 1
    n0 = 0.5
    rng = np.random.default_rng(42)
    x = rng.choice([-1, 1], size=200)
    # Forney model: y = conv(x, v) + white_noise
    from scipy.signal import lfilter
    signal = lfilter(v, [1.0], x)
    noise = rng.standard_normal(len(x)) * np.sqrt(n0 / 2.0)
    y = signal + noise
    result = cs_equalizer_bpsk(y, v, n0=n0, nu=2)
    assert np.all(np.isfinite(result.llr))
    assert np.mean(np.sign(result.llr) == x) > 0.5


def test_fd_shortening_params():
    """Frequency-domain shortening parameters should produce valid results."""
    v = np.array([0.8, 0.5, 0.2, 0.1])
    n0 = 0.5
    Z_w, g_diag, g_off, n_fft = compute_shortened_params_from_v(v, n0, nu=2)
    assert np.all(np.isfinite(Z_w))
    assert g_diag > 0
    assert len(g_off) == 2


def test_fd_ber_decreases_with_snr():
    """BER from FD-based shortened BCJR should decrease with increasing SNR."""
    v = np.array([0.8, 0.5, 0.2, 0.1])
    n_sig = 500
    n_fft = 1024
    bers = []
    for n0 in [1.0, 0.3, 0.05]:
        Z_w, g_diag, g_off, _ = compute_shortened_params_from_v(v, n0, nu=2, n_fft=n_fft)
        rng = np.random.default_rng(42)
        errors = 0
        total = 0
        from scipy.signal import lfilter
        for _ in range(10):
            x = rng.choice([-1, 1], size=n_sig)
            signal = lfilter(v, [1.0], x)
            noise = rng.standard_normal(len(x)) * np.sqrt(n0 / 2.0)
            y = signal + noise
            z = apply_shortening_filter_fft(y, Z_w, n_fft)
            result = shortened_bcjr_bpsk(z, g_diag, g_off, nu=2)
            errors += np.sum(np.sign(result.llr) != x)
            total += len(x)
        bers.append(errors / total)
    assert bers[-1] < bers[0]
