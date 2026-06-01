from __future__ import annotations

import numpy as np
import pytest

from ftn.pulse import compute_g, generate_rrc
from ftn.channel import ftn_awgn_channel, ftn_awgn_channel_complex
from ftn.modulation import bpsk_modulate, qam16_modulate
from ftn.equalizers.simplified_mbcjr import simplified_mbcjr_bpsk, simplified_mbcjr_qam16
from ftn.equalizers.ungerboeck_mbcjr import ungerboeck_mbcjr_bpsk


def _make_g():
    t, h = generate_rrc(beta=0.3, span=15, sps=128)
    return compute_g(t, h, tau=0.5, isi_len=3)


class TestSimplifiedMbcjrBpsk:
    def test_survivors_within_limit(self):
        g = _make_g()
        rng = np.random.default_rng(42)
        bits = rng.integers(0, 2, size=50, dtype=np.uint8)
        x = bpsk_modulate(bits)
        y = ftn_awgn_channel(x, g, n0=0.5, rng=rng)
        result = simplified_mbcjr_bpsk(y, g, n0=0.5, m_states=4, future_len=3)
        assert all(c <= 4 for c in result.survivor_counts)

    def test_llr_finite(self):
        g = _make_g()
        rng = np.random.default_rng(43)
        bits = rng.integers(0, 2, size=30, dtype=np.uint8)
        x = bpsk_modulate(bits)
        y = ftn_awgn_channel(x, g, n0=1.0, rng=rng)
        result = simplified_mbcjr_bpsk(y, g, n0=1.0, m_states=4, future_len=3)
        assert np.all(np.isfinite(result.bit_llr))
        assert result.bit_llr.shape == (30,)

    def test_close_to_original_with_large_m(self):
        g = _make_g()
        rng = np.random.default_rng(44)
        bits = rng.integers(0, 2, size=40, dtype=np.uint8)
        x = bpsk_modulate(bits)
        y = ftn_awgn_channel(x, g, n0=0.3, rng=rng)
        n_states = 2 ** 3
        orig = ungerboeck_mbcjr_bpsk(y, g, n0=0.3, m_states=n_states, future_len=5)
        simp = simplified_mbcjr_bpsk(y, g, n0=0.3, m_states=n_states, future_len=5)
        corr = np.corrcoef(orig.llr, simp.bit_llr)[0, 1]
        assert corr > 0.8, f"Corr={corr:.3f}, expected > 0.8"

    def test_ber_decreases_with_snr(self):
        g = _make_g()
        bers = []
        for snr_db in [1.0, 3.0, 5.0]:
            rng = np.random.default_rng(45 + int(snr_db * 10))
            bits = rng.integers(0, 2, size=200, dtype=np.uint8)
            x = bpsk_modulate(bits)
            n0 = 1.0 / (10.0 ** (snr_db / 10.0))
            y = ftn_awgn_channel(x, g, n0=n0, rng=rng)
            result = simplified_mbcjr_bpsk(y, g, n0=n0, m_states=4, future_len=5)
            hard = (result.bit_llr < 0).astype(np.uint8)
            bers.append(float(np.mean(hard != bits)))
        assert bers[-1] < bers[0], f"BER should decrease: {bers}"


class TestSimplifiedMbcjrQam16:
    def test_llr_shape(self):
        g = _make_g()
        rng = np.random.default_rng(60)
        bits = rng.integers(0, 2, size=20, dtype=np.uint8)
        x = qam16_modulate(bits)
        y = ftn_awgn_channel_complex(x, g, n0=0.5, rng=rng)
        result = simplified_mbcjr_qam16(y, g, n0=0.5, m_states=4, future_len=2)
        assert result.bit_llr.shape == (20,)
        assert np.all(np.isfinite(result.bit_llr))

    def test_survivors_within_limit(self):
        g = _make_g()
        rng = np.random.default_rng(61)
        bits = rng.integers(0, 2, size=16, dtype=np.uint8)
        x = qam16_modulate(bits)
        y = ftn_awgn_channel_complex(x, g, n0=0.5, rng=rng)
        result = simplified_mbcjr_qam16(y, g, n0=0.5, m_states=4, future_len=2)
        assert all(c <= 4 for c in result.survivor_counts)

    def test_bit_llr_sign_no_noise(self):
        g = _make_g()
        rng = np.random.default_rng(62)
        bits = rng.integers(0, 2, size=8, dtype=np.uint8)
        x = qam16_modulate(bits)
        y = ftn_awgn_channel_complex(x, g, n0=1e-6, rng=rng)
        result = simplified_mbcjr_qam16(y, g, n0=1e-6, m_states=16, future_len=2)
        hard = (result.bit_llr < 0).astype(np.uint8)
        # At very low noise, most bits should be correct
        ber = float(np.mean(hard != bits))
        assert ber < 0.5, f"Noise-free BER={ber}"
