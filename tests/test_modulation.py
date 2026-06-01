import numpy as np

from ftn.modulation import bpsk_hard_bits_from_llr, bpsk_modulate


def test_bpsk_modulate_handles_uint8_without_underflow():
    bits = np.array([0, 1, 1, 0], dtype=np.uint8)

    symbols = bpsk_modulate(bits)

    np.testing.assert_array_equal(symbols, np.array([1.0, -1.0, -1.0, 1.0]))


def test_bpsk_hard_bits_follow_plus_over_minus_llr_convention():
    llr = np.array([4.0, -0.1, 0.0])

    hard = bpsk_hard_bits_from_llr(llr)

    np.testing.assert_array_equal(hard, np.array([0, 1, 0], dtype=np.uint8))
