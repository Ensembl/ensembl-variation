#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: build_vr_memory_hints.sh [--trace TRACE_FILE] [--work-dir WORK_DIR] RUN_DIR [OUTPUT_FILE]

Build Variant Recoder memory hints from a completed MaveDB run.

Defaults:
  TRACE_FILE  RUN_DIR/reports/trace.txt
  WORK_DIR    RUN_DIR/work
  OUTPUT_FILE RUN_DIR/vr_memory_hints.tsv
EOF
}

trace_file=""
work_dir=""
positional=()

while (( $# )); do
  case "$1" in
    -h|--help)
      usage
      exit 0
      ;;
    --trace|--work-dir)
      if (( $# < 2 )); then
        echo "ERROR: $1 requires a path" >&2
        exit 2
      fi
      if [[ $1 == "--trace" ]]; then
        trace_file=$2
      else
        work_dir=$2
      fi
      shift 2
      ;;
    -*)
      echo "ERROR: unknown option: $1" >&2
      exit 2
      ;;
    *)
      positional+=("$1")
      shift
      ;;
  esac
done

if (( ${#positional[@]} < 1 || ${#positional[@]} > 2 )); then
  usage >&2
  exit 2
fi

run_dir=$(realpath "${positional[0]}")
trace_file=${trace_file:-"${run_dir}/reports/trace.txt"}
work_dir=${work_dir:-"${run_dir}/work"}
output_file=${positional[1]:-"${run_dir}/vr_memory_hints.tsv"}

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
