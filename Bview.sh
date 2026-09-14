#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"

case_name="$(basename "${script_dir}")"
velocity_name="$(basename "$(dirname "${script_dir}")")"
case_id="${velocity_name}_${case_name}"
declared_case_id="$(awk '/^#define[[:space:]]+CASE_ID[[:space:]]/ {gsub(/\"/, "", $3); print $3; exit}' constants.h)"
if [ "${declared_case_id}" != "${case_id}" ]; then
    echo "[ERROR] Folder case ID '${case_id}' does not match constants.h '${declared_case_id}'."
    exit 1
fi

viewer="BpostBview_${case_id}"

qcc -Wall -O2 -disable-dimensions BpostBview.c -o "${viewer}" -L$BASILISK/gl -lglutils -lfb_tiny -lm
"./${viewer}" tb0.00 ts0.01 te1.00
