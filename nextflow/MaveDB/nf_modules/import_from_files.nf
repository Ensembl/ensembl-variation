process import_from_files {
  tag { urn }
  env.MAVEDB_URN = { urn }
  env.STEP       = 'import_from_files'

  input:
    val urn

  output:
    tuple val(urn), file("mappings.json"), file("scores.csv")

  script:
  """
  #!/usr/bin/env bash
  set +e

  log() {
    local ts; ts="$(date -Is)"
    >&2 echo "[$ts][MaveDB][URN=${MAVEDB_URN:-na}][STEP=${STEP:-na}][REASON=$1][SUBID=${2:-na}] ${3:-}"
  }

  log "stage_start"

  import_from_files.sh ${urn} ${params.mappings_path} ${params.scores_path}
  rc=$?
  if [ "$rc" -ne 0 ]; then
    log "stager_nonzero_exit" "na" "rc=$rc"
  fi

  if [ ! -s mappings.json ]; then
    echo "{}" > mappings.json
    log "missing_mappings_after_staging; wrote empty JSON"
  fi

  if [ ! -s scores.csv ]; then
    echo "" > scores.csv
    log "missing_scores_after_staging; wrote empty CSV"
  fi

  log "stage_complete"
  """
}