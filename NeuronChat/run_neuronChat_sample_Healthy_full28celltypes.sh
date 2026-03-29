#!/bin/bash
# ══════════════════════════════════════════════════════════════════
# LAUNCHER: CellChat Analysis
# Usage: bash launch_cellchat_AD_Dyslexia.sh
#        nohup bash launch_cellchat_AD_Dyslexia.sh > run_AD_Dyslexia.log 2>&1 &
# ══════════════════════════════════════════════════════════════════

# ── Thread control (prevent oversubscription) ─────────────────────
# Without this, BLAS spawns threads per worker = server killer
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1   # Apple/Intel MKL fallback
export NUMEXPR_NUM_THREADS=1      # numpy fallback (just in case)

# ── Temp directory ────────────────────────────────────────────────
export TMPDIR="/mnt/nfs/CX000008_DS1/projects/btanasa/tmp"
mkdir -p "$TMPDIR"

# ── Memory limits ─────────────────────────────────────────────────
ulimit -v unlimited    # virtual memory
ulimit -s unlimited    # stack size

# ── Log system state at launch ────────────────────────────────────
echo "╔══════════════════════════════════════════╗"
echo "║  CELLCHAT LAUNCHER                       ║"
echo "╚══════════════════════════════════════════╝"
echo "Started     : $(date)"
echo "Host        : $(hostname)"
echo "Cores avail : $(nproc)"
echo "Total RAM   : $(free -h | awk '/^Mem/{print $2}')"
echo "Free RAM    : $(free -h | awk '/^Mem/{print $4}')"
echo "TMPDIR      : $TMPDIR"
echo "TMPDIR free : $(df -h $TMPDIR | awk 'NR==2{print $4}')"
echo ""

# ── Script path ───────────────────────────────────────────────────
R_SCRIPT="/mnt/nfs/CX000008_DS1/projects/jaeyeon/fastq_file_Dyslexia_r1/zanalysis_bogdan/samples.merged_AG_Harmony_res0.x_anno_031826/script_neuronChat_Healthy_28celltypes.R"

# script_cellChat_Healthy_full28celltypes.R

echo "R script    : $R_SCRIPT"
echo ""

# ── Run ───────────────────────────────────────────────────────────
Rscript "$R_SCRIPT"

# ── Log completion ────────────────────────────────────────────────
EXIT_CODE=$?
echo ""
echo "══════════════════════════════════════════"
if [ $EXIT_CODE -eq 0 ]; then
  echo "✓ COMPLETED SUCCESSFULLY"
else
  echo "❌ FAILED with exit code: $EXIT_CODE"
fi
echo "Ended       : $(date)"
echo "Free RAM    : $(free -h | awk '/^Mem/{print $4}')"
echo "══════════════════════════════════════════"

exit $EXIT_CODE
