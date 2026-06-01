import numpy as np

from ftn.coding.conv import conv_bcjr_decode, conv_encode_75


def test_conv_encode_75_matches_manual_short_sequence():
    bits = np.array([1, 0, 1, 1], dtype=np.uint8)

    encoded = conv_encode_75(bits)

    np.testing.assert_array_equal(
        encoded,
        np.array([1, 1, 1, 0, 0, 0, 0, 1], dtype=np.uint8),
    )


def test_conv_bcjr_decode_has_zero_ber_with_strong_channel_llrs():
    bits = np.array([1, 0, 1, 1, 0, 0], dtype=np.uint8)
    encoded = conv_encode_75(bits)
    code_llr = np.where(encoded == 0, 30.0, -30.0)

    result = conv_bcjr_decode(code_llr)
    hard = (result.info_llr < 0.0).astype(np.uint8)

    np.testing.assert_array_equal(hard, bits)
    assert np.all(np.isfinite(result.info_llr))
    assert np.all(np.isfinite(result.code_llr))
    assert result.code_llr.shape == code_llr.shape
