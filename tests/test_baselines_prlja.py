"""Tests for Paper [14] Prlja M-BCJR baseline algorithms."""

import numpy as np
import pytest

from ftn.baselines.forney_model import (
    forney_channel,
    min_phase_from_pulse,
    spectral_factorize,
)
from ftn.baselines.prlja_mbcjr import (
    PrljaMbcjrResult,
    prlja_backup_mbcjr_bpsk,
    prlja_mbcjr_bpsk,
)
from ftn.pulse import generate_rrc


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_test_signal(n: int = 100, seed: int = 42):
    """Create a simple BPSK test signal with Forney model."""
    rng = np.random.default_rng(seed)
    x = rng.choice([-1, 1], size=n).astype(np.float64)
    return x, rng


def _simple_g():
    """A simple ISI autocorrelation for testing."""
    return {0: 1.0, 1: 0.4, -1: 0.4, 2: 0.1, -2: 0.1}


# ---------------------------------------------------------------------------
# Simple M-BCJR tests
# ---------------------------------------------------------------------------

class TestPrljaSimpleMbcjr:
    """Tests for prlja_mbcjr_bpsk (Simple Detection M-BCJR)."""

    def test_finite_llrs(self):
        """All LLRs should be finite."""
        x, rng = _make_test_signal(50)
        g = _simple_g()
        v = spectral_factorize(g)
        y = forney_channel(x, v, n0=0.5, rng=rng)
        result = prlja_mbcjr_bpsk(y, v, n0=0.5, M=4)
        assert np.all(np.isfinite(result.llr))

    def test_output_length(self):
        """LLR length should match input length."""
        x, rng = _make_test_signal(30)
        g = _simple_g()
        v = spectral_factorize(g)
        y = forney_channel(x, v, n0=0.5, rng=rng)
        result = prlja_mbcjr_bpsk(y, v, n0=0.5, M=4)
        assert len(result.llr) == 30
        assert len(result.survivor_counts) == 30

    def test_survivor_count_le_M(self):
        """Survivor count at each position should not exceed M."""
        x, rng = _make_test_signal(50)
        g = _simple_g()
        v = spectral_factorize(g)
        y = forney_channel(x, v, n0=0.5, rng=rng)
        M = 4
        result = prlja_mbcjr_bpsk(y, v, n0=0.5, M=M)
        for sc in result.survivor_counts:
            assert sc <= M

    def test_return_type(self):
        """Should return PrljaMbcjrResult."""
        x, rng = _make_test_signal(20)
        g = _simple_g()
        v = spectral_factorize(g)
        y = forney_channel(x, v, n0=0.5, rng=rng)
        result = prlja_mbcjr_bpsk(y, v, n0=0.5, M=4)
        assert isinstance(result, PrljaMbcjrResult)

    def test_ber_decreases_with_snr(self):
        """BER should decrease as SNR increases (N0 decreases)."""
        n_bits = 200
        n0_high = 2.0
        n0_low = 0.1

        g = _simple_g()
        v = spectral_factorize(g)
        x, rng = _make_test_signal(n_bits, seed=42)

        y_high = forney_channel(x, v, n0=n0_high, rng=np.random.default_rng(42))
        y_low = forney_channel(x, v, n0=n0_low, rng=np.random.default_rng(42))

        res_high = prlja_mbcjr_bpsk(y_high, v, n0=n0_high, M=8)
        res_low = prlja_mbcjr_bpsk(y_low, v, n0=n0_low, M=8)

        dec_high = (res_high.llr >= 0).astype(float) * 2 - 1
        dec_low = (res_low.llr >= 0).astype(float) * 2 - 1

        ber_high = np.mean(dec_high != x)
        ber_low = np.mean(dec_low != x)

        assert ber_low <= ber_high

    def test_m_increase_no_degradation(self):
        """Larger M should not produce worse BER."""
        n_bits = 100
        g = _simple_g()
        v = spectral_factorize(g)
        x, rng = _make_test_signal(n_bits, seed=99)
        y = forney_channel(x, v, n0=0.5, rng=rng)

        res_m4 = prlja_mbcjr_bpsk(y, v, n0=0.5, M=4)
        res_m16 = prlja_mbcjr_bpsk(y, v, n0=0.5, M=16)

        dec_m4 = (res_m4.llr >= 0).astype(float) * 2 - 1
        dec_m16 = (res_m16.llr >= 0).astype(float) * 2 - 1

        ber_m4 = np.mean(dec_m4 != x)
        ber_m16 = np.mean(dec_m16 != x)

        assert ber_m16 <= ber_m4 + 0.05  # small tolerance for stochastic effects

    def test_no_isi_channel(self):
        """With v=[1], should work like a simple AWGN channel."""
        x, rng = _make_test_signal(50)
        v = np.array([1.0])
        n0 = 0.5
        y = forney_channel(x, v, n0=n0, rng=rng)
        result = prlja_mbcjr_bpsk(y, v, n0=n0, M=4)
        # BER should be reasonable for this SNR
        dec = (result.llr >= 0).astype(float) * 2 - 1
        ber = np.mean(dec != x)
        assert ber < 0.2  # very loose, just checking it works

    def test_with_prior(self):
        """Should accept a-priori LLRs."""
        x, rng = _make_test_signal(30)
        g = _simple_g()
        v = spectral_factorize(g)
        y = forney_channel(x, v, n0=0.5, rng=rng)
        la = np.zeros(30, dtype=np.float64)
        result = prlja_mbcjr_bpsk(y, v, n0=0.5, M=4, la=la)
        assert np.all(np.isfinite(result.llr))

    def test_invalid_M(self):
        """M <= 0 should raise ValueError."""
        x = np.array([1.0, -1.0])
        v = np.array([1.0])
        with pytest.raises(ValueError):
            prlja_mbcjr_bpsk(x, v, n0=0.5, M=0)


# ---------------------------------------------------------------------------
# Backup M-BCJR tests
# ---------------------------------------------------------------------------

class TestPrljaBackupMbcjr:
    """Tests for prlja_backup_mbcjr_bpsk."""

    def test_finite_llrs(self):
        """Backup M-BCJR should produce finite LLRs."""
        x, rng = _make_test_signal(50)
        g = _simple_g()
        v = spectral_factorize(g)
        y = forney_channel(x, v, n0=0.5, rng=rng)
        result = prlja_backup_mbcjr_bpsk(y, v, n0=0.5, M=4, M_B=2)
        assert np.all(np.isfinite(result.llr))

    def test_output_length(self):
        """Output length should match input length."""
        x, rng = _make_test_signal(30)
        g = _simple_g()
        v = spectral_factorize(g)
        y = forney_channel(x, v, n0=0.5, rng=rng)
        result = prlja_backup_mbcjr_bpsk(y, v, n0=0.5, M=4, M_B=2)
        assert len(result.llr) == 30

    def test_smooth_flag(self):
        """Smooth flag should not crash and produce finite LLRs."""
        x, rng = _make_test_signal(50)
        g = _simple_g()
        v = spectral_factorize(g)
        y = forney_channel(x, v, n0=0.5, rng=rng)
        result = prlja_backup_mbcjr_bpsk(y, v, n0=0.5, M=4, M_B=2, smooth=True)
        assert np.all(np.isfinite(result.llr))

    def test_backup_improves_reserve_values(self):
        """Backup should replace reserve values with computed LLRs."""
        x, rng = _make_test_signal(100)
        g = {0: 1.0, 1: 0.5, -1: 0.5, 2: 0.3, -2: 0.3, 3: 0.1, -3: 0.1}
        v = spectral_factorize(g)
        y = forney_channel(x, v, n0=0.8, rng=rng)

        simple = prlja_mbcjr_bpsk(y, v, n0=0.8, M=2)
        backup = prlja_backup_mbcjr_bpsk(y, v, n0=0.8, M=2, M_B=4)

        # Check that reserve values (|llr|~5.0) are handled
        LAMBDA = 5.0
        n_reserve_simple = np.sum(np.abs(np.abs(simple.llr) - LAMBDA) < 0.01)
        n_reserve_backup = np.sum(np.abs(np.abs(backup.llr) - LAMBDA) < 0.01)

        # Backup should not have more reserves than simple
        assert n_reserve_backup <= n_reserve_simple

    def test_real_ftn_minphase_llrs_are_finite(self):
        """Real FTN minimum-phase front-end should not leak infinities."""
        rng = np.random.default_rng(42)
        x = rng.choice([-1.0, 1.0], size=200)
        t, h = generate_rrc(beta=0.3, span=15, sps=128)
        v = min_phase_from_pulse(t, h, tau=0.5, g0=1.0)
        y = forney_channel(x, v, n0=0.6325, rng=rng)

        result = prlja_backup_mbcjr_bpsk(y, v, n0=0.6325, M=4, M_B=2, smooth=True)

        assert np.all(np.isfinite(result.llr))
