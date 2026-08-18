#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: build_vr_memory_hints.sh RUN_DIR [OUTPUT_FILE]

Build Variant Recoder memory hints from a completed MaveDB run.

RUN_DIR must contain:
  trace.txt
  work/

OUTPUT_FILE defaults to RUN_DIR/vr_memory_hints.tsv.
EOF
}

if [[ ${1:-} == "-h" || ${1:-} == "--help" ]]; then
  usage
  exit 0
fi

if (( $# < 1 || $# > 2 )); then
  usage >&2
  exit 2
fi

run_dir=$(realpath "$1")
trace_file="${run_dir}/trace.txt"
work_dir="${run_dir}/work"
output_file=${2:-"${run_dir}/vr_memory_hints.tsv"}

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
generator="${script_dir}/../bin/build_vr_memory_hints.py"

if [[ ! -f "${trace_file}" ]]; then
  echo "ERROR: trace file not found: ${trace_file}" >&2
  exit 1
fi

if [[ ! -d "${work_dir}" ]]; then
  echo "ERROR: work directory not found: ${work_dir}" >&2
  exit 1
fi

if [[ ! -x "${generator}" ]]; then
  echo "ERROR: memory-hint generator is not executable: ${generator}" >&2
  exit 1
fi

output_dir=$(dirname "${output_file}")
if [[ ! -d "${output_dir}" ]]; then
  echo "ERROR: output directory not found: ${output_dir}" >&2
  exit 1
fi

temporary_output="${output_file}.tmp.$$"
trap 'rm -f "${temporary_output}"' EXIT

python3 "${generator}" \
  --trace "${trace_file}" \
  --work-dir "${work_dir}" \
  --output "${temporary_output}"

mv "${temporary_output}" "${output_file}"
trap - EXIT

echo "Memory hints written to: ${output_file}"
echo "Use on the next run with: --vr_memory_hints ${output_file}"
