#!/usr/bin/env nextflow

process crossmap {
  tag "$id"
  container 'quay.io/biocontainers/crossmap:0.7.0--pyhdfd78af_0'
  memory '4GB'
  time '1h'

  input:
    path vcf
    each id
    path chain_dir
    path fasta_dir

  output:
    tuple val(id), path("*.vcf"), path("*.unmap"), emit: vcf
    path("*_report.txt"), emit: report

  script:
    def chain = "${chain_dir}/*${id}*"
    def fasta = "${fasta_dir}/*${id}*.bgz"
    def out   = "${vcf.simpleName}_${id}.vcf"
  """
  CrossMap vcf ${chain} ${vcf} ${fasta} ${out} 2> ${id}_report.txt
  """
}

process tabix {
  tag "$id"
  memory '4GB'
  time '1h'

  publishDir "${params.out_dir}", mode: 'move'

  input:
    tuple val(id), path(vcf), path(unmap)
  output:
    path("*.vcf.gz*")

  """
  for file in $vcf $unmap; do
    sort -k1,1d -k2,2n \${file} -o \${file}
    bgzip \${file}
    tabix -p vcf \${file}.gz
  done
  """
}

process report {
  memory '1GB'
  time '1h'

  publishDir "${params.out_dir}", mode: 'move'

  input:
    path report
  output:
    path "${params.report}"

  """
  echo "#ID\tFailed\tTotal\tPercentage" > ${params.report}
  for i in `echo ${report} | tr " " "\n" | sort`; do
    total=`grep -Eo "Total entries: .*" \$i | grep -Eo "[0-9]+"`
    failed=`grep -Eo "Failed to map: .*" \$i | grep -Eo "[0-9]+"`
    percentage=`bc <<< "scale=2; \$failed * 100 / \$total"`
    echo "\${i/_report.txt/}\t\$failed\t\$total\t\$percentage" >> ${params.report}
  done
  """
}
