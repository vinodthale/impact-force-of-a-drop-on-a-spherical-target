#!/usr/bin/env bash
# Set number of MPI processes (change here as needed)
NUM_PROCS=6

case_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
case_name="$(basename "${case_dir}")"
velocity_name="$(basename "$(dirname "${case_dir}")")"
case_id="${velocity_name}_${case_name}"
cd "${case_dir}"

declared_case_id="$(awk '/^#define[[:space:]]+CASE_ID[[:space:]]/ {gsub(/\"/, "", $3); print $3; exit}' constants.h)"
if [ "${declared_case_id}" != "${case_id}" ]; then
    echo "[ERROR] Folder case ID '${case_id}' does not match constants.h '${declared_case_id}'."
    exit 1
fi

executable="Bdropimpact_${case_id}"

# Define log file
LOGFILE="simulation_${case_id}_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOGFILE") 2>&1

echo "===== Simulation started at $(date) ====="
echo "Case ID: $case_id"
echo "Using $NUM_PROCS MPI processes"
echo "Working directory: $(pwd)"
echo "Log file: $LOGFILE"
echo

# Step 1: Generate MPI-enabled source with disabled dimensions
echo "[INFO] Generating source code with qcc..."
if ! qcc -source -D_MPI=1 -disable-dimensions Bdropimpact.c; then
    echo "[ERROR] qcc failed to generate source code."
    exit 1
fi
echo "[OK] qcc generated source _Bdropimpact.c"
echo

# Step 2: Compile
echo "[INFO] Compiling with mpicc..."
if ! mpicc -O2 -Wall -std=c99 -D_MPI=1 -D_GNU_SOURCE=1 -D_FORTIFY_SOURCE=0 _Bdropimpact.c -o "$executable" -lm; then
    echo "[ERROR] Compilation failed."
    exit 1
fi
echo "[OK] Compilation succeeded."
echo

# Step 3: Run the simulation
echo "[INFO] Running simulation with $NUM_PROCS MPI processes..."
START_TIME=$(date +%s)

if ! mpirun -np "$NUM_PROCS" "./${executable}"; then
    echo "[ERROR] Simulation failed during execution."
    exit 1
fi

END_TIME=$(date +%s)
echo
echo "===== Simulation finished at $(date) ====="
echo "Total time: $((END_TIME - START_TIME)) seconds"







