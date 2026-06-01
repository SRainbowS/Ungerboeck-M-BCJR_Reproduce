import numpy as np

from ftn.channel import build_isi_matrix, ftn_filter_output, sample_colored_noise


def test_filter_output_matches_toeplitz_matrix_product():
    g = {0: 1.0, 1: 0.25, -1: 0.25, 2: -0.08, -2: -0.08}
    x = np.array([1.0, -1.0, 1.0, 1.0, -1.0])

    G = build_isi_matrix(g, len(x))
    y_matrix = G @ x
    y_filter = ftn_filter_output(x, g)

    np.testing.assert_allclose(y_filter, y_matrix, atol=1e-12)


def test_colored_noise_empirical_covariance_matches_ungerboeck_model():
    rng = np.random.default_rng(1234)
    g = {0: 1.0, 1: 0.18, -1: 0.18, 2: 0.05, -2: 0.05}
    n0 = 0.7
    n = 4

    samples = sample_colored_noise(g, n=n, n0=n0, rng=rng, size=30000)
    empirical = np.cov(samples, rowvar=False, bias=True)
    expected = (n0 / 2.0) * build_isi_matrix(g, n)

    np.testing.assert_allclose(empirical, expected, atol=1.8e-2)
