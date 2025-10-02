#!/bin/bash

log() {
  local ts; ts="$(date -Is)"
  >&2 echo "[$ts][MaveDB][URN=${MAVEDB_URN:-${urn:-na}}][STEP=${STEP:-import_from_files}][REASON=$1][SUBID=${2:-na}] ${3:-}"
}

log "import_start" "na" "msg=import_from_files-processing-URN:'${1:-na}'"

if [[ $# < 3 ]]
then 
      log "bad_args" "na" "got=$# expected>=3"
      exit 1
fi
urn=$1
mappings_path=$2
scores_path=$3

# make URN available for downstream uniform logs if not already set
export MAVEDB_URN="${MAVEDB_URN:-$urn}"
export STEP="${STEP:-import_from_files}"

# Locate the mapping file using the original urn (with colons)
mapping_file=$(find ${mappings_path} -type f -iname "*${urn}*.json" | head -n 1)

# Replace colons with hyphens for searching the scores directory
score_urn=$(echo "${urn}" | sed 's/:/-/g')
score_file=$(find ${scores_path} -type f -iname "*${score_urn}.scores.csv" | head -n 1)

# Check if the mapping and scores files exist in the user-provided directories
if [ ! -f "${mapping_file}" ]; then
  log "mapping_missing" "na" "searched=${mappings_path} pattern=*${urn}*.json"
  exit 1
fi
if [ ! -f "${score_file}" ]; then
  log "scores_missing" "na" "searched=${scores_path} pattern=*${score_urn}.scores.csv"
  exit 1
fi

# If the files exist, report this
log "mapping_found" "na" "src=${mapping_file}"
log "scores_found" "na" "src=${score_file} sanitized=${score_urn}"

# Copy mappings file to the working dir
cp "${mapping_file}" mappings.json

# Check if the file contains any lines starting with "tmp:"
if grep -q '^tmp:' ${score_file}; then

  log "scores_tmp_ids_detected"

  # Get the file's base name (e.g., "urn-mavedb-00000001-a-1")
  prefix=$(basename "${score_file}")
  prefix=${prefix%.scores.csv}
  log "scores_prefix" "na" "prefix=${prefix}"

  # Process the file with awk:
  # 1. Substitute "urn-mavedb-" with "urn:mavedb:"
  # 2. Print the columns
  # 3. For each subsequent row, substitute the first field: replace the pattern starting with "tmp:" up to and including the "#" with file prefix (file name) and "#"
  awk -F, -v OFS=, -v prefix="$prefix" '
    BEGIN {
      gsub(/^urn-mavedb-/, "urn:mavedb:", prefix)
    }
    NR==1 { print }
    NR>1 { sub(/^tmp:[^#]+#/, prefix "#" , $1); print }
  ' "${score_file}" >"scores.csv"

else
  # If the file does not contain any lines starting with "tmp:", copy the file as is
  log "scores_ids_ok_copy"
  cp "$score_file" ./scores.csv
fi

# If contents of pwd is mappings.json and scores.csv, then the files are copied successfully - check
if [ -f mappings.json ] && [ -f scores.csv ]; then
  log "stage_ok"
  log "stage_listing" "na" "pwd=$(pwd)"
  ls -l mappings.json scores.csv 1>&2 || true
else
  log "stage_failed"
  log "stage_listing" "na" "pwd=$(pwd)"
  ls -l 1>&2 || true
  exit 1
fi