from __future__ import annotations

import numpy as np
import pytest

from ftn.coding.turbo import (
    rsc_encode_4state,
    rsc_encode_8state,
    turbo_encode,
    turbo_decode,
    s_random_interleaver,
    TurboEncodeResult,
)


class TestRsc4State:
    def test_known_output(self):
        info = np.array([1, 0, 1, 1, 0], dtype=np.uint8)
        sys, parity, fin = rsc_encode_4state(info)
        assert np.array_equal(sys, info)
        assert parity.shape == info.shape
        assert 0 <= fin <= 3

    def test_zero_input_terminates(self):
        info = np.array([0, 0, 0, 0], dtype=np.uint8)
        sys, parity, fin = rsc_encode_4state(info)
        assert fin == 0


class TestRsc8State:
    def test_known_output(self):
        info = np.array([1, 0, 1, 1, 0], dtype=np.uint8)
        sys, parity, fin = rsc_encode_8state(info)
        assert np.array_equal(sys, info)
        assert parity.shape == info.shape
        assert 0 <= fin <= 7

    def test_zero_input_terminates(self):
        info = np.array([0, 0, 0, 0, 0], dtype=np.uint8)
        sys, parity, fin = rsc_encode_8state(info)
        assert fin == 0


class TestTurboEncode:
    def test_output_length(self):
        K = 100
        info = np.zeros(K, dtype=np.uint8)
        pi = s_random_interleaver(K, seed=0)
        result = turbo_encode(info, pi)
        # Layout: [sys(K+5), p1(K+5), p2(K+5)] where 5=mu1+mu2
        assert result.coded_bits.size == 3 * (K + 5)
        assert result.n_info == K
        assert result.n_tail1 == 2
        assert result.n_tail2 == 3

    def test_all_zero(self):
        K = 50
        info = np.zeros(K, dtype=np.uint8)
        pi = s_random_interleaver(K, seed=0)
        result = turbo_encode(info, pi)
        assert result.coded_bits.dtype == np.uint8


class TestSRandomInterleaver:
    def test_valid_permutation(self):
        K = 200
        pi = s_random_interleaver(K, seed=42)
        assert pi.shape == (K,)
        assert np.array_equal(np.sort(pi), np.arange(K))

    def test_deterministic(self):
        K = 100
        pi1 = s_random_interleaver(K, seed=7)
        pi2 = s_random_interleaver(K, seed=7)
        assert np.array_equal(pi1, pi2)


class TestTurboDecode:
    def test_noiseless_roundtrip(self):
        K = 100
        rng = np.random.default_rng(123)
        info = rng.integers(0, 2, size=K, dtype=np.uint8)
        pi = s_random_interleaver(K, seed=0)
        enc = turbo_encode(info, pi)

        # BPSK on all coded bits: 0->+1, 1->-1
        symbols = 1.0 - 2.0 * enc.coded_bits.astype(float)
        llr = 20.0 * symbols

        # Layout: [sys(K+5), p1(K+5), p2(K+5)]
        n = K + 5
        sys_llr = llr[:n]
        p1_llr = llr[n:2 * n]
        p2_llr = llr[2 * n:3 * n]

        result = turbo_decode(sys_llr, p1_llr, p2_llr, pi, max_iterations=20)
        hard = (result.info_llr < 0).astype(np.uint8)
        ber = np.mean(hard != info)
        assert ber == 0.0, f"BER={ber}, expected 0"
        assert result.converged

    def test_returns_result(self):
        K = 50
        pi = s_random_interleaver(K, seed=1)
        enc = turbo_encode(np.zeros(K, dtype=np.uint8), pi)
        symbols = 1.0 - 2.0 * enc.coded_bits.astype(float)
        llr = 5.0 * symbols
        n = K + 5
        sys_llr = llr[:n]
        p1_llr = llr[n:2 * n]
        p2_llr = llr[2 * n:3 * n]
        result = turbo_decode(sys_llr, p1_llr, p2_llr, pi, max_iterations=5)
        assert result.info_llr.shape == (K,)
        assert result.n_iterations >= 1
