process map_scores_to_HGVSp_variants {
  tag { urn }
  env.MAVEDB_URN = { urn }
  env.STEP       = 'map_scores_hgvsp'

  input:  tuple val(urn), path(mappings), path(scores), path(metadata), path(vr)
  output: tuple val(urn), path('map_*.tsv')

  memory { mappings.size() * 4.B + 1.GB }

  script:
  def round = params.round ? "--round ${params.round}" : ""
  def script_name = params.from_files ? "map_scores_to_variants_fromfiles.py" : "map_scores_to_variants.py"

  """
  #!/usr/bin/env bash
  set +e

  log() { local ts; ts="\$(date -Is)"; >&2 echo "[\$ts][MaveDB][URN=\${MAVEDB_URN:-na}][STEP=\${STEP:-na}][REASON=\$1][SUBID=\${2:-na}] \${3:-}"; }

  log "stage_start"

  ${script_name} --urn ${urn} \\
                 --scores ${scores} \\
                 --mappings ${mappings} \\
                 --metadata ${metadata} \\
                 --vr ${vr} \\
                 ${round} \\
                 --output map_${urn}.tsv
  rc=\$?
  if [ "\$rc" -ne 0 ]; then
    log "mapper_nonzero_exit" "na" "rc=\$rc"
  fi

  if [ ! -s map_${urn}.tsv ]; then
      : > map_${urn}.tsv
      log "mapper_empty_output_created"
  fi

  log "stage_complete"
  """
}

process map_scores_to_HGVSg_variants {
  tag { urn }
  env.MAVEDB_URN = { urn }
  env.STEP       = 'map_scores_hgvs_g'

  input:  tuple val(urn), path(mappings), path(scores), path(metadata), val(hgvs)
  output: tuple val(urn), path(metadata), path('*map_*.tsv')

  memory { mappings.size() * 2.B + 1.GB }

  script:
  def round = params.round ? "--round ${params.round}" : ""
  def script_name = params.from_files ? "map_scores_to_variants_fromfiles.py" : "map_scores_to_variants.py"

  """
  #!/usr/bin/env bash
  set +e

  log() { local ts; ts="\$(date -Is)"; >&2 echo "[\$ts][MaveDB][URN=\${MAVEDB_URN:-na}][STEP=\${STEP:-na}][REASON=\$1][SUBID=\${2:-na}] \${3:-}"; }

  log "stage_start"

  ${script_name} --urn ${urn} \\
                 --scores ${scores} \\
                 --mappings ${mappings} \\
                 --metadata ${metadata} \\
                 ${round} \\
                 --output map_${urn}.tsv
  rc=\$?
  if [ "\$rc" -ne 0 ]; then
    log "mapper_nonzero_exit" "na" "rc=\$rc"
  fi

  if [ ! -s map_${urn}.tsv ]; then
      : > map_${urn}.tsv
      log "mapper_empty_output_created"
  fi

  log "stage_complete"
  """
}