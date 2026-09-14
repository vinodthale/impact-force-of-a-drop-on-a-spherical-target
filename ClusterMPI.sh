#!/bin/bash
#SBATCH -J V5.00_0_deg
#SBATCH -N 1
#SBATCH --ntasks-per-node=64
#SBATCH -p tyhcnormal
#SBATCH -t 7-00:00:00
#SBATCH -o sim-%j.out
#SBATCH -e sim-%j.err
#SBATCH --exclusive

# ==========================================
# Basilisk MPI SLURM Run Script
# ==========================================

set -eo pipefail

# ==========================================
# Load environment
# ==========================================
module purge

source /work/home/acs5n14dk2/soft/basilisk-26.08.05-gcc12.2.0/scripts/env.sh
module load compiler/openmpi/4.1.8-gcc12.2.0

# ==========================================
# Logging
# ==========================================
LOGFILE="runlog_${SLURM_JOB_ID}.log"
exec > >(tee -a "$LOGFILE") 2>&1

# ==========================================
# Diagnostics
# ==========================================
echo "=========================================="
echo "Basilisk MPI Job"
echo "=========================================="
echo "SLURM Job ID : $SLURM_JOB_ID"
echo "Node List    : $SLURM_JOB_NODELIST"
echo "Task Count   : $SLURM_NTASKS"
echo "Start Time   : $(date)"
echo "Working Dir  : $PWD"
echo "Executable   : ./Bdropimpact"
echo "=========================================="
echo

echo "[INFO] MPI environment"
echo "qcc   : $(command -v qcc)"
echo "mpicc : $(command -v mpicc)"
echo "srun  : $(command -v srun)"
echo

# ==========================================
# Pre-run validation
# ==========================================
if [ ! -f "./Bdropimpact" ]; then
    echo "[ERROR] Executable ./Bdropimpact not found."
    echo "Compile first:"
    echo "bash compile_Basilisk.sh"
    exit 1
fi

# ==========================================
# Run simulation
# ==========================================
echo "[INFO] Starting simulation with $SLURM_NTASKS MPI ranks..."

START_TIME=$(date +%s)

srun --mpi=pmix_v3 ./Bdropimpact

END_TIME=$(date +%s)

# ==========================================
# Completion
# ==========================================
echo
echo "[OK] Simulation completed successfully"
echo "Runtime: $((END_TIME - START_TIME)) seconds"
echo "Finished at: $(date)"
