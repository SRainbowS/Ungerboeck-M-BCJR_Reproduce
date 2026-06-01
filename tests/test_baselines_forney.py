"""Tests for Forney model spectral factorization and channel."""

import numpy as np
import pytest

from ftn.baselines.forney_model import (
    forney_branch_metric,
    forney_channel,
    spectral_factorize,
)


class TestSpectralFactorization:
    """Tests for spectral_factorize."""

    def test_unit_g(self):
        """g = {0: 1.0} should produce v = [1.0]."""
        g = {0: 1.0}
        v = spectral_factorize(g)
        assert len(v) == 1
        np.testing.assert_allclose(v, [1.0], atol=1e-12)

    def test_energy_conservation(self):
        """sum(v^2) should equal g[0]."""
        g = {0: 1.0, 1: 0.3, 2: 0.1}
        v = spectral_factorize(g)
        np.testing.assert_allclose(np.sum(v ** 2), g[0], atol=1e-10)

    def test_min_phase_property(self):
        """v[0] should be the largest absolute value in v."""
        g = {0: 1.0, 1: 0.5, 2: 0.2}
        v = spectral_factorize(g)
        assert abs(v[0]) >= np.max(np.abs(v)) - 1e-12

    def test_reconstruction_autocorrelation(self):
        """Autocorrelation of v should reproduce g."""
        g = {0: 1.0, 1: 0.4, 2: 0.15}
        v = spectral_factorize(g)
        L = len(v) - 1
        for lag in range(L + 1):
            auto = sum(v[k] * v[k + lag] for k in range(len(v) - lag))
            np.testing.assert_allclose(auto, g[lag], atol=1e-10)

    def test_symmetric_g(self):
        """g with symmetric lags (as autocorrelation) should work."""
        g = {0: 1.3, 1: 0.6, -1: 0.6, 2: 0.2, -2: 0.2}
        v = spectral_factorize(g)
        assert len(v) == 3  # L = max lag = 2, so len = 3
        # Verify energy
        np.testing.assert_allclose(np.sum(v ** 2), g[0], atol=1e-10)

    def test_longer_isi(self):
        """Longer ISI sequence should factorize correctly."""
        g = {0: 1.0, 1: 0.5, 2: 0.3, 3: 0.1, -1: 0.5, -2: 0.3, -3: 0.1}
        v = spectral_factorize(g)
        L = len(v) - 1
        assert L == 3
        for lag in range(L + 1):
            auto = sum(v[k] * v[k + lag] for k in range(len(v) - lag))
            np.testing.assert_allclose(auto, g[lag], atol=1e-10)


class TestForneyChannel:
    """Tests for forney_channel."""

    def test_unit_channel(self):
        """v=[1] should produce y=x (no noise)."""
        rng = np.random.default_rng(42)
        x = rng.choice([-1, 1], size=20).astype(float)
        y = forney_channel(x, np.array([1.0]), n0=0.0)
        np.testing.assert_allclose(y, x, atol=1e-12)

    def test_noiseless_is_channel(self):
        """Noiseless channel output should match causal convolution."""
        rng = np.random.default_rng(123)
        x = rng.choice([-1, 1], size=50).astype(float)
        v = np.array([1.0, 0.5, 0.2])
        y = forney_channel(x, v, n0=0.0)
        # Verify against manual convolution
        from scipy.signal import lfilter
        expected = lfilter(v, [1.0], x)
        np.testing.assert_allclose(y, expected, atol=1e-12)

    def test_noisy_output_shape(self):
        """Output shape should match input shape."""
        rng = np.random.default_rng(7)
        x = rng.choice([-1, 1], size=30).astype(float)
        v = np.array([1.0, 0.3])
        y = forney_channel(x, v, n0=0.5, rng=rng)
        assert y.shape == x.shape

    def test_noisy_variance(self):
        """Noise variance should be approximately N0/2."""
        rng = np.random.default_rng(99)
        x = np.ones(10000)
        v = np.array([1.0])
        n0 = 1.0
        y = forney_channel(x, v, n0=n0, rng=rng)
        noise = y - x
        # noise variance should be ~ n0/2 = 0.5
        np.testing.assert_allclose(np.var(noise), n0 / 2.0, atol=0.05)


class TestForneyBranchMetric:
    """Tests for forney_branch_metric."""

    def test_zero_residual(self):
        """Perfect match should give -(0)^2/N0 = 0 plus prior."""
        bm = forney_branch_metric(y_n=1.0, label=1.0, n0=0.5)
        np.testing.assert_allclose(bm, 0.0, atol=1e-12)

    def test_with_prior(self):
        """Prior should be added."""
        bm = forney_branch_metric(y_n=1.0, label=1.0, n0=0.5, log_prior=2.0)
        np.testing.assert_allclose(bm, 2.0, atol=1e-12)

    def test_negative_metric_for_mismatch(self):
        """Mismatched y and label should give negative metric."""
        bm = forney_branch_metric(y_n=1.0, label=-1.0, n0=0.5)
        assert bm < 0.0
        expected = -(2.0 ** 2) / 0.5
        np.testing.assert_allclose(bm, expected, atol=1e-12)
