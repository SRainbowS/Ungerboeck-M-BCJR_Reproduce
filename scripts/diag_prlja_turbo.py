"""Diagnostic: why does prlja turbo give BER ~0.49?

Run a single frame through the turbo loop and inspect intermediate values.
"""
import sys
import numpy as np
from pathlib import Path

PROJECT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT))

from ftn.baselines.forney_model import forney_channel, min_phase_from_pulse
from ftn.baselines.prlja_mbcjr import prlja_mbcjr_bpsk, prlja_backup_mbcjr_bpsk
from ftn.channel import ftn_awgn_channel
from ftn.coding.conv import conv_bcjr_decode, conv_encode_75
from ftn.equalizers.ungerboeck_mbcjr import ungerboeck_mbcjr_bpsk
from ftn.modulation import bpsk_modulate
from ftn.pulse import compute_g, generate_rrc

tau = 0.5
isi_len = 7
snr_db = 5.0
n_bits = 1000
seed = 42

# Setup channel
t, h = generate_rrc(beta=0.3, span=15, sps=128)
g = compute_g(t, h, tau=tau, isi_len=isi_len)
v = min_phase_from_pulse(t, h, tau=tau, g0=g[0])

print(f"v length: {len(v)}, v[:5] = {v[:5]}")
print(f"v energy: {np.sum(v**2):.4f}, g[0] = {g[0]:.4f}")

code_rate = 0.5
n0_es = 1.0 / (10.0 ** ((snr_db + 10*np.log10(code_rate)) / 10.0))
print(f"n0_es = {n0_es:.4f}")

rng = np.random.default_rng(seed)

# Generate data
info_bits = rng.integers(0, 2, size=n_bits, dtype=np.uint8)
encoded = conv_encode_75(info_bits)
print(f"encoded length: {len(encoded)}, info_bits length: {len(info_bits)}")

interleaver = np.random.default_rng(seed + 1000).permutation(2 * n_bits)
inv_perm = np.argsort(interleaver)

int_encoded = encoded[interleaver]
symbols = bpsk_modulate(int_encoded)

print(f"symbols: min={symbols.min():.2f}, max={symbols.max():.2f}")

# Forney channel
y_forney = forney_channel(symbols, v, n0=n0_es, rng=rng)
# Ungerboeck channel
y_ung = ftn_awgn_channel(symbols, g, n0=n0_es, rng=rng)

print(f"y_forney: mean={y_forney.mean():.4f}, std={y_forney.std():.4f}")
print(f"y_ung:    mean={y_ung.mean():.4f}, std={y_ung.std():.4f}")

# --- Test 1: Uncoded detection (no turbo loop) ---
print("\n=== Test 1: Uncoded detection ===")

# prlja
det_prlja = prlja_backup_mbcjr_bpsk(y_forney, v, n0=n0_es, M=4, M_B=2, smooth=True)
decisions_prlja = (det_prlja.llr >= 0).astype(int)  # LLR>0 -> bit 0 -> symbol +1
bits_prlja = ((1 - decisions_prlja) > 0).astype(np.uint8)  # bit 1 for symbol -1
ber_prlja_uncoded = np.mean(bits_prlja != int_encoded)
print(f"prlja M=4 uncoded BER: {ber_prlja_uncoded:.4f}")
print(f"  LLR stats: mean={det_prlja.llr.mean():.4f}, std={det_prlja.llr.std():.4f}, "
      f"min={det_prlja.llr.min():.2f}, max={det_prlja.llr.max():.2f}")

# ungerboeck
initial_state = tuple(1.0 for _ in range(isi_len))
det_ung = ungerboeck_mbcjr_bpsk(y_ung, g, n0=n0_es, isi_len=isi_len,
                                  m_states=4, future_len=3, initial_state=initial_state)
decisions_ung = (det_ung.llr >= 0).astype(int)
bits_ung = ((1 - decisions_ung) > 0).astype(np.uint8)
ber_ung_uncoded = np.mean(bits_ung != int_encoded)
print(f"ung M=4 uncoded BER: {ber_ung_uncoded:.4f}")

# --- Test 2: Turbo loop ---
print("\n=== Test 2: Turbo equalization loop ===")

llr_clip = 20.0
turbo_iters = 5

for algo_name, detector_fn in [
    ("prlja_M4", lambda y, prior: prlja_backup_mbcjr_bpsk(
        y, v, n0=n0_es, M=4, M_B=2, smooth=True, la=prior)),
    ("ung_M4", lambda y, prior: ungerboeck_mbcjr_bpsk(
        y, g, n0=n0_es, la=prior, isi_len=isi_len,
        m_states=4, future_len=3, initial_state=initial_state)),
]:
    print(f"\n--- {algo_name} ---")
    y = y_forney if "prlja" in algo_name else y_ung
    detector_prior = np.zeros(y.size, dtype=float)

    for it in range(turbo_iters + 1):
        det = detector_fn(y, detector_prior)
        det_llr = np.clip(det.llr, -llr_clip, llr_clip)
        det_ext = np.clip(det_llr - detector_prior, -llr_clip, llr_clip)

        # Stats on detector extrinsic
        print(f"  Iter {it}: det_llr mean={det_llr.mean():.4f} std={det_llr.std():.4f} "
              f"| det_ext mean={det_ext.mean():.4f} std={det_ext.std():.4f} "
              f"| prior mean={detector_prior.mean():.4f} std={detector_prior.std():.4f}")

        # De-interleave
        det_ext_deint = det_ext[inv_perm]

        # Decode
        dec = conv_bcjr_decode(det_ext_deint)
        info_llr = np.clip(dec.info_llr, -llr_clip, llr_clip)
        code_llr = np.clip(dec.code_llr, -llr_clip, llr_clip)

        # Decoder extrinsic
        dec_ext = np.clip(code_llr - det_ext_deint, -llr_clip, llr_clip)
        print(f"         info_llr mean={info_llr.mean():.4f} std={info_llr.std():.4f} "
              f"| code_llr mean={code_llr.mean():.4f} std={code_llr.std():.4f} "
              f"| dec_ext mean={dec_ext.mean():.4f} std={dec_ext.std():.4f}")

        # Final decisions (from info_llr)
        hard_info = (info_llr < 0).astype(np.uint8)  # LLR<0 -> bit 1
        ber = np.mean(hard_info != info_bits)
        print(f"         BER after iter {it}: {ber:.4f} ({np.sum(hard_info != info_bits)}/{len(info_bits)})")

        # Interleave decoder extrinsic for next iteration
        detector_prior = dec_ext[interleaver]

print("\nDone.")
