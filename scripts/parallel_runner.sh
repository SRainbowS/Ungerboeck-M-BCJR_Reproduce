#!/usr/bin/env bash
# ============================================================================
# FTN M-BCJR 并行蒙特卡洛仿真启动脚本
#
# 用法:
#   bash scripts/parallel_runner.sh fig7           # Fig.7, 默认参数
#   bash scripts/parallel_runner.sh fig7 52        # Fig.7, 52 并行进程
#   bash scripts/parallel_runner.sh fig8 16        # Fig.8, 16 并行进程
#   bash scripts/parallel_runner.sh fig7 4 3000 20 # 自定义 K=3000, 20帧/SNR
#
# 每个 SNR 点分配给一个独立进程，结果写入 results/ 子目录，
# 最后由 merge_results.py 合并。
# ============================================================================
set -euo pipefail

FIGURE="${1:-fig7}"
N_PROC="${2:-52}"
K="${3:-6000}"
FRAMES="${4:-100}"

# 限制每个 NumPy 进程只用 1 线程，避免 OpenBLAS 线程过度抢占
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

PYTHON="python"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_DIR"

echo "============================================"
echo "FTN M-BCJR Parallel Monte Carlo Simulation"
echo "============================================"
echo "Figure:    $FIGURE"
echo "Processes: $N_PROC"
echo "K (bits):  $K"
echo "Frames:    $FRAMES per SNR point"
echo "============================================"

if [ "$FIGURE" = "fig7" ]; then
    SNR_RANGE=$(python -c "import numpy as np; print(' '.join(f'{x:.2f}' for x in np.arange(2.0, 6.01, 0.25)))")
    SCRIPT="scripts/reproduce_fig7.py"
    OUTDIR="results/fig7_tau_05"
    TAU="0.5"
    ISI_LEN=7
    TURBO_ITERS=5
    CONFIGS="M2_L3:m2_l3:2:3 M2_L5:m2_l5:2:5"
elif [ "$FIGURE" = "fig8" ]; then
    SNR_RANGE=$(python -c "import numpy as np; print(' '.join(f'{x:.2f}' for x in np.arange(2.5, 6.51, 0.25)))")
    SCRIPT="scripts/reproduce_fig8.py"
    OUTDIR="results/fig8_tau_035"
    TAU="0.35"
    ISI_LEN=10
    TURBO_ITERS=15
    CONFIGS="M8_L5:m8_l5:8:5 M8_L7:m8_l7:8:7"
else
    echo "Unknown figure: $FIGURE (use fig7 or fig8)"
    exit 1
fi

SNR_LIST=($SNR_RANGE)
N_SNR=${#SNR_LIST[@]}
echo "SNR points: $N_SNR"
echo "Output dir: $OUTDIR"
echo ""

# 创建临时目录存放各进程的部分结果
PARTIAL_DIR="$OUTDIR/partial_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$PARTIAL_DIR"

PIDS=()
TASK_ID=0
N_RUNNING=0

# 为每个 SNR 点和配置启动一个独立进程
for CONFIG_STR in $CONFIGS; do
    IFS=':' read -r LABEL DIRNAME M_STATES FUTURE_LEN <<< "$CONFIG_STR"

    for SNR_DB in "${SNR_LIST[@]}"; do
        # 节流：用计数器+wait -n 替代失效的 jobs -r 子 shell 方案
        while [ $N_RUNNING -ge $N_PROC ]; do
            wait -n 2>/dev/null || true
            N_RUNNING=$((N_RUNNING - 1))
        done

        SNR_SAFE=$(echo "$SNR_DB" | tr '.' '_')
        PARTIAL_OUT="$PARTIAL_DIR/${LABEL}_snr_${SNR_SAFE}.csv"
        SEED=$((20260507 + TASK_ID))

        echo "[Task $TASK_ID] $LABEL @ Eb/N0=${SNR_DB} dB → $PARTIAL_OUT"

        $PYTHON -c "
import sys, csv, numpy as np
sys.path.insert(0, '$PROJECT_DIR')
from ftn.pulse import generate_rrc, compute_g
from ftn.channel import ftn_awgn_channel
from ftn.coding.conv import conv_encode_75, conv_bcjr_decode
from ftn.equalizers.ungerboeck_mbcjr import ungerboeck_mbcjr_bpsk
from ftn.modulation import bpsk_modulate

TAU = $TAU
ROLLOFF = 0.3
PULSE_SPAN = 15
SPS = 128
ISI_LEN = $ISI_LEN
M_STATES = $M_STATES
FUTURE_LEN = $FUTURE_LEN
TURBO_ITERS = $TURBO_ITERS
LLR_CLIP = 20.0
CODE_RATE = 0.5
K = $K
EB_N0 = $SNR_DB
ES_N0 = EB_N0 + 10 * np.log10(CODE_RATE)
N0 = 1.0 / (10 ** (ES_N0 / 10.0))
SEED = $SEED
FRAMES = $FRAMES
MIN_ERRORS = 200

t, h = generate_rrc(beta=ROLLOFF, span=PULSE_SPAN, sps=SPS)
g = compute_g(t, h, tau=TAU, isi_len=ISI_LEN)
interleaver = np.random.default_rng(SEED + 1000).permutation(2 * K)
inv_perm = np.argsort(interleaver)
initial_state = tuple(1.0 for _ in range(ISI_LEN))

rng = np.random.default_rng(SEED)
bit_errors = 0
total_bits = 0
for fi in range(FRAMES):
    info_bits = rng.integers(0, 2, size=K, dtype=np.uint8)
    encoded = conv_encode_75(info_bits)
    int_encoded = encoded[interleaver]
    symbols = bpsk_modulate(int_encoded)
    y = ftn_awgn_channel(symbols, g, n0=N0, rng=rng)
    detector_prior = np.zeros(y.size, dtype=float)
    for _ in range(TURBO_ITERS + 1):
        det = ungerboeck_mbcjr_bpsk(y, g, n0=N0, la=detector_prior,
            isi_len=ISI_LEN, m_states=M_STATES, future_len=FUTURE_LEN,
            initial_state=initial_state)
        det_llr = np.clip(det.llr, -LLR_CLIP, LLR_CLIP)
        det_ext = np.clip(det_llr - detector_prior, -LLR_CLIP, LLR_CLIP)
        det_ext_deint = det_ext[inv_perm]
        dec = conv_bcjr_decode(det_ext_deint)
        info_llr = np.clip(dec.info_llr, -LLR_CLIP, LLR_CLIP)
        code_llr = np.clip(dec.code_llr, -LLR_CLIP, LLR_CLIP)
        dec_ext = np.clip(code_llr - det_ext_deint, -LLR_CLIP, LLR_CLIP)
        detector_prior = dec_ext[interleaver]
    hard_info = (info_llr < 0).astype(np.uint8)
    bit_errors += int(np.sum(hard_info != info_bits))
    total_bits += int(info_bits.size)
    if bit_errors >= MIN_ERRORS and fi >= 5:
        break

ber = bit_errors / total_bits if total_bits > 0 else 0.0
with open('$PARTIAL_OUT', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['eb_n0_db', 'label', 'm_states', 'future_len', 'bit_errors', 'total_bits', 'ber'])
    w.writerow([$SNR_DB, '$LABEL', $M_STATES, $FUTURE_LEN, bit_errors, total_bits, f'{ber:.10g}'])
print(f'[$LABEL] Eb/N0={EB_N0:.2f} dB  BER={ber:.6g}  ({bit_errors}/{total_bits})')
" &

        PIDS+=($!)
        TASK_ID=$((TASK_ID + 1))
        N_RUNNING=$((N_RUNNING + 1))
    done
done

echo ""
echo "Waiting for all $TASK_ID tasks to finish..."
for PID in "${PIDS[@]}"; do
    wait "$PID" 2>/dev/null || echo "Warning: PID $PID exited with error"
done

echo ""
echo "All tasks completed. Merging results..."
$PYTHON scripts/merge_results.py "$PARTIAL_DIR" "$OUTDIR" "$FIGURE"

echo ""
echo "============================================"
echo "Done! Results in $OUTDIR"
echo "============================================"
