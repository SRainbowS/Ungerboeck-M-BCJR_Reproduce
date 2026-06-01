import numpy as np

from ftn.channel import ftn_filter_output
from ftn.equalizers.full_bcjr import brute_force_bpsk_llr, full_bcjr_bpsk


def test_full_bcjr_matches_brute_force_map_for_short_sequence():
    g = {0: 1.0, 1: 0.35, -1: 0.35, 2: -0.12, -2: -0.12}
    x_true = np.array([1.0, -1.0, -1.0, 1.0, -1.0])
    y = ftn_filter_output(x_true, g) + np.array([0.05, -0.1, 0.02, 0.07, -0.03])
    la = np.array([0.2, -0.3, 0.1, 0.0, 0.4])
    n0 = 0.8

    result = full_bcjr_bpsk(y, g, n0=n0, la=la, isi_len=2, initial_state=(1.0, 1.0))
    brute = brute_force_bpsk_llr(y, g, n0=n0, la=la, isi_len=2, initial_state=(1.0, 1.0))

    np.testing.assert_allclose(result.llr, brute, atol=1e-10)
