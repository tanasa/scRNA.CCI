#!/bin/bash
#
# run_scDiffCom.sh  —  Launch scDiffCom R script with proper environment
#

# ── Thread control (prevent oversubscription) ───────────────
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLAS_NUM_THREADS=1

# ── Temp directory ───────────────────────────────────────────
export TMPDIR="/mnt/nfs/CX000008_DS1/projects/btanasa/tmp"
mkdir -p "$TMPDIR"

# ── Memory limits ────────────────────────────────────────────
ulimit -v unlimited
ulimit -s unlimited

# ── Log system state at launch ───────────────────────────────
echo "Started     : $(date)"
echo "Host        : $(hostname)"
echo "Cores avail : $(nproc)"
echo "Free RAM    : $(free -h | awk '/^Mem/{print $4}')"
echo "TMPDIR      : $TMPDIR"
echo ""

# ── Run R script ─────────────────────────────────────────────
Rscript /mnt/nfs/script_run_scDiffCom_comparison.R

echo "Finished    : $(date)"
