#!/usr/bin/env nextflow

process copy_to_rapid_release_ftp {
  memory '1GB'
  time '1h'

  errorStrategy 'terminate'
  cache false

  input:
    tuple val(id), path(vcf)
    val lookup

  output:
    stdout

  script:
    def assembly = !params.keep_id && lookup[id] ? lookup[id] : id
    def outdir = params.rr_root + "/" + params.rr_path
  """
  assembly=${assembly}
  for file in ${vcf}; do
    become ensrapid rsync -av `realpath \${file}` ${outdir}
  done

  # Add or update checksums
  become ensrapid rm -f ${outdir}/md5sum* ${outdir}/CHECKSUMS
  become ensrapid bash -c "cd ${outdir} && md5sum * > CHECKSUMS"
  """
}
