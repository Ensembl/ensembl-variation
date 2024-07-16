#!/usr/bin/env nextflow

process copy_to_rapid_release_ftp {
  memory '1GB'
  time '10m'
  
  errorStrategy 'terminate'

  input:
    tuple val(id), path(vcf)
    val lookup

  output:
    stdout

  script:
    def assembly = !params.keep_id && lookup[id] ? lookup[id] : id
  """
  assembly=${assembly}
  for file in ${vcf}; do
    become ensrapid rsync -av `realpath \${file}` ${params.rr_root}/${params.rr_path}
  done
  """
}
