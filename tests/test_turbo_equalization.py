import numpy as np

from ftn.channel import ftn_awgn_channel
from ftn.coding.conv import conv_encode_75
from ftn.modulation import bpsk_modulate
from ftn.pulse import compute_g, generate_rrc
from ftn.turbo_equalization import turbo_equalize_conv_bpsk


def test_turbo_equalization_keeps_llrs_finite_and_does_not_worsen_high_snr_ber():
    info_bits = np.array([1, 0, 1, 1, 0, 0, 1, 0], dtype=np.uint8)
    encoded = conv_encode_75(info_bits)
    symbols = bpsk_modulate(encoded)
    t, h = generate_rrc(beta=0.3, span=15, sps=128)
    g = compute_g(t, h, tau=0.8, isi_len=2)
    rng = np.random.default_rng(9)
    y = ftn_awgn_channel(symbols, g, n0=0.05, rng=rng)

    result = turbo_equalize_conv_bpsk(
        y,
        g,
        n0=0.05,
        info_bits=info_bits,
        turbo_iters=2,
        isi_len=2,
        m_states=4,
        future_len=3,
        llr_clip=20.0,
    )

    bers = [item.info_ber for item in result.iterations]
    assert bers[-1] <= bers[0]
    assert np.all(np.isfinite(result.final_info_llr))
    assert np.all(np.isfinite(result.final_detector_llr))
