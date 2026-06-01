from __future__ import annotations

from types import SimpleNamespace

import numpy as np

import scripts.reproduce_fig11 as fig11


def test_fig11_ftn_frame_interleaves_code_bits_and_deinterleaves_llrs(monkeypatch):
    info_bits = np.array([1, 0, 1, 1], dtype=np.uint8)
    encoded = np.array([0, 0, 1, 1, 0, 1, 1, 0], dtype=np.uint8)
    interleaver = np.array([2, 5, 0, 7, 1, 6, 3, 4])
    inv_perm = np.argsort(interleaver)

    detector_ext = np.array([3.0, -4.0, 5.0, -6.0, 7.0, -8.0, 9.0, -10.0])
    decoder_code_llr = detector_ext[inv_perm] + np.array(
        [0.5, -0.5, 1.5, -1.5, 2.5, -2.5, 3.5, -3.5]
    )
    expected_detector_prior = decoder_code_llr - detector_ext[inv_perm]

    detector_calls = []
    decoder_calls = []

    def fake_encode(bits):
        assert np.array_equal(bits, info_bits)
        return encoded

    def fake_modulate(bits):
        assert np.array_equal(bits, encoded[interleaver])
        return np.zeros(2, dtype=complex)

    def fake_channel(symbols, g, n0, rng):
        assert symbols.shape == (2,)
        return np.zeros(2, dtype=complex)

    def fake_detector(y, g, n0, la, isi_len, m_states, future_len):
        detector_calls.append(la.copy())
        if len(detector_calls) == 1:
            assert np.array_equal(la, np.zeros(8))
        else:
            assert np.allclose(la, expected_detector_prior[interleaver])
        return SimpleNamespace(bit_llr=la + detector_ext)

    def fake_decode(code_llr):
        decoder_calls.append(code_llr.copy())
        assert np.allclose(code_llr, detector_ext[inv_perm])
        return SimpleNamespace(
            info_llr=np.array([-1.0, 1.0, -1.0, -1.0]),
            code_llr=decoder_code_llr,
        )

    monkeypatch.setattr(fig11, "conv_encode_75", fake_encode)
    monkeypatch.setattr(fig11, "qam16_modulate", fake_modulate)
    monkeypatch.setattr(fig11, "ftn_awgn_channel_complex", fake_channel)
    monkeypatch.setattr(fig11, "simplified_mbcjr_qam16", fake_detector)
    monkeypatch.setattr(fig11, "conv_bcjr_decode", fake_decode)

    errors, total = fig11._ftn_turbo_frame_16qam(
        info_bits,
        {0: 1.0},
        n0=1.0,
        isi_len=3,
        m_states=4,
        future_len=3,
        turbo_iters=1,
        llr_clip=20.0,
        interleaver=interleaver,
        rng=np.random.default_rng(0),
    )

    assert (errors, total) == (0, info_bits.size)
    assert len(detector_calls) == 2
    assert len(decoder_calls) == 2
