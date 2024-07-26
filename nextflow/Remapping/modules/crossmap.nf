#!/usr/bin/env nextflow

process crossmap {
  tag "$id"
  container 'quay.io/biocontainers/crossmap:0.7.0--pyhdfd78af_0'
  memory '4GB'
  time '1h'

  input:
    path vcf
    each id
    path chain_dir, stageAs: 'chain'
    path fasta_dir, stageAs: 'fasta'

  output:
    tuple val(id), path("*.vcf"), path("*.unmap"), emit: vcf
    path("*_report.txt"), emit: report

  script:
    def chain = "${chain_dir}/*${id}*.chain{,.gz,.bgz}"
    def fasta = "${fasta_dir}/*${id}*.{fa,fasta}{,.gz,.bgz}"
    def out   = "${vcf.simpleName}_${id}.vcf"
  """
  shopt -s nullglob # do not expand glob on non-matching files
  CrossMap vcf ${chain} ${vcf} ${fasta} ${out} 2> ${id}_report.txt
  """
}

process tabix {
  tag "$id"
  memory '4GB'
  time '1h'

  publishDir "${params.out_dir}", mode: 'copy'

  input:
    tuple val(id), path(vcf), path(unmap)
    val lookup

  output:
    tuple val(id), path("*.vcf.gz*")

  script:
    def assembly = !params.keep_id && lookup[id] ? lookup[id] : id
  """
  for file in $vcf $unmap; do
    final=\${file/$id/$assembly}
    (grep "^#" \${file} && grep -v "^#" \${file} | sort -k1,1d -k2,2n) > \${final}
    bgzip \${final}
    tabix -p vcf \${final}.gz
  done
  """
}

process report {
  memory '1GB'
  time '1h'

  publishDir "${params.out_dir}", mode: 'copy'

  input:
    path report
    val lookup
  output:
    path "${params.report}"

  script:
    def mapping = lookup.collect { "['${it.key}']=\"${it.value}\"" }.join(' ')
  """
  declare -A arr=(${mapping})
  for i in `echo ${report} | tr " " "\n" | sort`; do
    id=\${i/_report.txt/}
    total=`grep -Eo "Total entries: .*" \$i | grep -Eo "[0-9]+"`
    failed=`grep -Eo "Failed to map: .*" \$i | grep -Eo "[0-9]+"`
    percentage=`bc <<< "scale=2; \$failed * 100 / \$total"`

    # set assembly to empty string if no match is found
    assembly=`[ -v arr[\$id] ] && echo \${arr[\$id]} || echo ''`
    echo "\$id\t\$assembly\t\$failed\t\$total\t\$percentage" >> ${params.report}
  done

  # sort by percentage of failed variants (descending)
  sort -k5,5nr ${params.report} -o ${params.report}

  # add header
  sed -ie '1i\\#ID\tAssembly\tFailed\tTotal\tPercentage' ${params.report}
  """
}
