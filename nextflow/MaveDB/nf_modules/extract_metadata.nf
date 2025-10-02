// extractMetadata parses and grabs the URN ID's entry from the pre-downloaded maveDB main.json file,
// which contains metadata for all URN ID's.
// It then reformats it to match the JSON format expected by the rest of the pipeline which 
// was originally written to handle a different JSON structure. It uses extract_metadata.py to do this. 
process extract_metadata {
  tag { urn }
  env.MAVEDB_URN = { urn }
  env.STEP       = 'extract_metadata'

  input:
    val urn
    path metadata_file

  output:
    tuple val(urn), file("metadata.json"), file("LICENCE.txt")

  script:
  """
  #!/usr/bin/env bash
  set +e

  log() {
    local ts; ts="\$(date -Is)"
    >&2 echo "[\$ts][MaveDB][URN=\${MAVEDB_URN:-na}][STEP=\${STEP:-na}][REASON=\$1][SUBID=\${2:-na}] \${3:-}"
  }

  log "stage_start"

  extract_metadata.py --metadata_file ${metadata_file} --urn "${urn}"
  rc=\$?
  if [ "\$rc" -ne 0 ]; then
    log "extractor_nonzero_exit" "na" "rc=\$rc"
  fi

  # If metadata.json is missing or empty, create fallback files (unchanged functionality)
  if [ ! -s metadata.json ]; then
      echo "{}" > metadata.json
      echo ""  > LICENCE.txt
      log "missing_metadata_after_extractor; wrote empty JSON/LICENCE"
  fi

  log "stage_complete"
  """
}