#!/usr/bin/env bash
# ============================================================================
# 服务器部署与运行指南
#
# 1. 将以下文件/目录复制到服务器:
#    ftn/           (核心算法包)
#    scripts/       (仿真脚本)
#    tests/         (可选，验证用)
#
# 2. 安装依赖:
#    pip install numpy scipy
#
# 3. 运行仿真 (见下方命令)
# ============================================================================
#
# === Fig.7: tau=0.5, M=2, turbo_iters=5 ===
#
# 全参数 (K=6000, 100帧/SNR, ~8-17小时 @ 8进程):
#   bash scripts/parallel_runner.sh fig7 8
#
# 缩减帧数 (K=6000, 30帧/SNR, ~5-10小时 @ 8进程):
#   bash scripts/parallel_runner.sh fig7 8 6000 30
#
# 缩减K (K=3000, 100帧/SNR, ~4-8小时 @ 8进程):
#   bash scripts/parallel_runner.sh fig7 8 3000 100
#
# === Fig.8: tau=0.35, M=8, turbo_iters=15 ===
#
# 注意: M=8 + future_len=7 极慢 (每帧约 1400s @ K=6000)
# 强烈建议先跑 Fig.7 确认正确性，再跑 Fig.8
#
# 推荐: K=3000, 50帧/SNR, 16进程 (~60小时):
#   bash scripts/parallel_runner.sh fig8 16 3000 50
#
# 最小验证: K=1000, 20帧/SNR, 16进程 (~10小时):
#   bash scripts/parallel_runner.sh fig8 16 1000 20
#
# === 合并结果 ===
# parallel_runner.sh 完成后会自动调用 merge_results.py
# 也可以手动合并:
#   python scripts/merge_results.py <partial_dir> <output_dir> fig7
#
# === 验证 ===
# 快速冒烟测试 (K=500, 3帧, ~3分钟):
#   python scripts/reproduce_fig7.py --quick
#   python scripts/reproduce_fig8.py --quick
#
# 运行单元测试:
#   python -m pytest tests/ -v
# ============================================================================
