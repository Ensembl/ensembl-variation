// extractMetadata parses and grabs the URN ID's entry from the pre-downloaded maveDB main.json file,
// which contains metadata for all URN ID's.
// It then reformats it to match the JSON format expected by the rest of the pipeline which 
// was originally written to handle a different JSON structure. It uses extract_metadata.py to do this. 
process extract_metadata {
  tag { urn }
  // publishDir "${params.output}", mode: 'copy', overwrite: true

  input:
    val urn
    path metadata_file

  output:
    tuple val(urn), file("metadata.json"), file("LICENCE.txt")

  script:
  """
  #!/usr/bin/env bash
  set +e
  export MAVEDB_URN='${urn}'
  export STEP='extract_metadata'

  log() {
    local ts; ts="\$(date -Is)"
    >&2 echo "[\$ts][MaveDB][URN=\${MAVEDB_URN:-na}][STEP=\${STEP:-na}][REASON=\$1][SUBID=\${2:-na}] \${3:-}"
  }

  log "stage_start" "na" "metadata_file=${metadata_file}"

  extract_metadata.py --metadata_file ${metadata_file} --urn "${urn}"
  rc=\$?
  if [ "\$rc" -ne 0 ]; then
    log "extractor_nonzero_exit" "na" "rc=\$rc"
  fi

  # If metadata.json is missing or empty, create fallback files
  if [ ! -s metadata.json ]; then
      echo "{}" > metadata.json
      echo ""  > LICENCE.txt
      log "missing_metadata_after_extractor; wrote empty JSON/LICENCE"
  fi

  meta_bytes=\$(wc -c < metadata.json 2>/dev/null || echo 0)
  lic_bytes=\$(wc -c < LICENCE.txt 2>/dev/null || echo 0)
  log "stage_complete" "na" "meta_bytes=\${meta_bytes} licence_bytes=\${lic_bytes} out=metadata.json,LICENCE.txt"
  """
}
