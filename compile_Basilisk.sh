#!/bin/bash

# ==========================================
# Basilisk MPI Compilation Script
# ==========================================

set -eo pipefail

echo "========================================="
echo "Basilisk MPI Compilation"
echo "========================================="

WORK_DIR=$(pwd)
echo "Working directory: $WORK_DIR"
echo

# ==========================================
# Load cluster environment
# ==========================================
echo "[INFO] Loading Basilisk and MPI environment..."

module purge

source /work/home/acs5n14dk2/soft/basilisk-26.08.05-gcc12.2.0/scripts/env.sh

module load compiler/openmpi/4.1.8-gcc12.2.0

echo "[OK] Environment loaded"
echo

# ==========================================
# Diagnostics
# ==========================================
echo "[INFO] Diagnostics"
echo "qcc   : $(command -v qcc)"
echo "mpicc : $(command -v mpicc)"
echo

qcc --version
echo
mpicc --version
echo

# ==========================================
# Clean previous build
# ==========================================
echo "[INFO] Cleaning old build files..."
rm -f _Bdropimpact.c Bdropimpact
echo

# ==========================================
# Generate Basilisk MPI source
# ==========================================
echo "[INFO] Generating MPI-compatible source..."

if ! qcc -source -D_MPI=1 -disable-dimensions Bdropimpact.c; then
    echo "[ERROR] qcc source generation failed."
    exit 1
fi

if [ ! -f "_Bdropimpact.c" ]; then
    echo "[ERROR] Generated source file missing."
    exit 1
fi

echo "[OK] Generated _Bdropimpact.c"
echo

# ==========================================
# Compile
# ==========================================
echo "[INFO] Compiling with MPI compiler..."

if ! mpicc -O3 -Wall -std=c99 -D_MPI=1 -D_GNU_SOURCE \
    _Bdropimpact.c \
    -o Bdropimpact \
    -lm; then
    echo "[ERROR] Compilation failed."
    exit 1
fi

if [ ! -f "Bdropimpact" ]; then
    echo "[ERROR] Executable not created."
    exit 1
fi
if [ ! -x "Bdropimpact" ]; then
    echo "[ERROR] Executable Bdropimpact is not executable."
    exit 1
fi

echo
echo "[OK] Compilation successful"
ls -lh Bdropimpact
echo
echo "Executable ready: $WORK_DIR/Bdropimpact"
