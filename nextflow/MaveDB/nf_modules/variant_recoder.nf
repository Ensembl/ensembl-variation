process run_variant_recoder {
  // Run Variant Recoder on a file with HGVS identifiers
  label 'bigmem'

  input:
    tuple val(urn), path(mappings), path(scores), path(metadata), path(hgvs), val(hint_lines), val(hint_peak_rss_gb)
  output:
    tuple val(urn), path(mappings), path(scores), path(metadata), path('vr.json')

  tag "${urn}"
  memory { 
    def n = file(hgvs.target).countLines()

    def fallback =
      n <= 100   ? 2.GB   :
      n <= 500   ? 8.GB   :
      n <= 1000  ? 24.GB  :
      n <= 2000  ? 48.GB  :
      n <= 5000  ? 96.GB  :
      n <= 10000 ? 180.GB :
      n <= 20000 ? 220.GB :
                   300.GB

    def base = fallback
    if (hint_lines > 0 && hint_peak_rss_gb > 0) {
      def scaledPeakGb = hint_peak_rss_gb * n / hint_lines
      def hintedGb = (long) Math.ceil((scaledPeakGb * 1.3 + 4) / 4) * 4
      base = [[hintedGb.GB, 8.GB].max(), 300.GB].min()
    }

    if (task.attempt == 1) {
      return base
    }

    if (!(task.exitStatus in [130, 137, 140])) {
      return base
    }

    def escalated =
      task.attempt == 2 ? [base * 2, fallback].max() :
      task.attempt == 3 ? [base * 4, fallback * 2, 200.GB].max() :
                          500.GB

    [escalated, 500.GB].min()
  }

  script:
  def bin = "${params.ensembl}/ensembl-vep"
  def reg = params.registry ? "--registry ${params.registry}" : ""
  def memorySource = hint_lines > 0 ? "historical" : "fallback"
  """
  #!/usr/bin/env bash
  set +e
  export MAVEDB_URN='${urn}'
  export STEP='variant_recoder'
  log() { local ts; ts=\$(date -Is); >&2 echo "[\$ts][MaveDB][URN=\${MAVEDB_URN:-na}][STEP=\${STEP:-na}][REASON=\$1][SUBID=\${2:-na}] \${3:-}"; }

  HGVS_LINES=\$(wc -l < "${hgvs}" 2>/dev/null || echo 0)
  log "vr_start" "na" "attempt=${task.attempt} requested_mem=${task.memory} hgvs_lines=\${HGVS_LINES} memory_source=${memorySource} hint_lines=${hint_lines} hint_peak_rss_gb=${hint_peak_rss_gb}"

  perl ${bin}/variant_recoder -i ${hgvs} --vcf_string ${reg} > vr.json 2> vr.stderr
  rc=\$?

  # always log final size of the output file
  BYTES=\$(wc -c < vr.json 2>/dev/null || echo 0)
  HSIZE=\$(du -h vr.json 2>/dev/null | cut -f1 || echo 0)
  log "vr_output_size" "na" "bytes=\${BYTES} human=\${HSIZE} file=vr.json"

  if [ "\$rc" -ne 0 ]; then
    log "vr_exit_nonzero" "na" "rc=\$rc attempt=${task.attempt}"
    exit "\$rc"
  fi

  if [ ! -s vr.json ]; then
    log "vr_empty_output" "na" "attempt=${task.attempt}"
    exit 1
  fi

  log "vr_exit_ok" "na" "rc=0"
  """
}
