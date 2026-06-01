import numpy as np

from ftn.pulse import compute_g, generate_rrc


def test_rrc_autocorrelation_is_normalized_and_symmetric():
    t, h = generate_rrc(beta=0.3, span=15, sps=128)

    energy = np.trapz(np.abs(h) ** 2, t)
    assert np.isclose(energy, 1.0, atol=1e-4)

    g = compute_g(t, h, tau=0.5, isi_len=6)
    assert np.isclose(g[0], 1.0, atol=1e-4)
    for lag in range(1, 7):
        assert np.isclose(g[lag], g[-lag], atol=1e-8)


def test_nyquist_spacing_has_small_remote_isi_taps():
    t, h = generate_rrc(beta=0.3, span=15, sps=128)
    g = compute_g(t, h, tau=1.0, isi_len=5)

    assert abs(g[0] - 1.0) < 1e-4
    assert max(abs(g[lag]) for lag in range(2, 6)) < 2e-3
