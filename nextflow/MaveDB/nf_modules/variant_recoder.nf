process run_variant_recoder {
  // Run Variant Recoder on a file with HGVS identifiers
  label 'bigmem'

  input:
    tuple val(urn), path(mappings), path(scores), path(metadata), path(hgvs)
  output:
    tuple val(urn), path(mappings), path(scores), path(metadata), path('vr.json')

  tag "${urn}"
  memory { 
    def want = file(hgvs.target).countLines() * 200.MB + 50.GB 
    def cap  = 200.GB
    [want, cap].min()
  }

  script:
  def bin = "${params.ensembl}/ensembl-vep"
  def reg = params.registry ? "--registry ${params.registry}" : ""
  """
  #!/usr/bin/env bash
  set +e
  export MAVEDB_URN='${urn}'
  export STEP='variant_recoder'
  log() { local ts; ts=\$(date -Is); >&2 echo "[\$ts][MaveDB][URN=\${MAVEDB_URN:-na}][STEP=\${STEP:-na}][REASON=\$1][SUBID=\${2:-na}] \${3:-}"; }

  HGVS_LINES=\$(wc -l < "${hgvs}" 2>/dev/null || echo 0)
  log "vr_start" "na" "attempt=${task.attempt} requested_mem=${task.memory} hgvs_lines=\${HGVS_LINES}"

  perl ${bin}/variant_recoder -i ${hgvs} --vcf_string ${reg} > vr.json 2> vr.stderr
  rc=\$?

  # always log final size of the output file
  BYTES=\$(wc -c < vr.json 2>/dev/null || echo 0)
  HSIZE=\$(du -h vr.json 2>/dev/null | cut -f1 || echo 0)
  log "vr_output_size" "na" "bytes=\${BYTES} human=\${HSIZE} file=vr.json"

  if [ "\$rc" -ne 0 ]; then
    log "vr_exit_nonzero" "na" "rc=\$rc attempt=${task.attempt}"
  else
    log "vr_exit_ok" "na" "rc=0"
  fi
  """
}
