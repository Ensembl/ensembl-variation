#!/usr/bin/env nextflow

process bcftools_liftover {

    tag "$id"
    container 'fairbrot/bcftools-liftover:latest'
    cpus    1
    memory  '4 GB'
    time    '1h'

    input:
        path vcf
        each id
        val  chainFile  // params.chain
        val  srcFa      // params.src_fasta
        val  dstFa      // params.dst_fasta

    output:
        tuple val(id), path("${vcf.simpleName}_${id}.vcf"), path("${vcf.simpleName}_${id}.unmap"), emit: vcf
        path ("${id}_report.txt"), emit: report

    script:
      def outVcf = "${vcf.simpleName}_${id}.vcf"
      def unmVcf = "${vcf.simpleName}_${id}.unmap"
    """
    set -euo pipefail

    # ---------- liftover ----------
    bcftools +liftover ${vcf} -Ov -o ${outVcf} -- \\
            -s ${srcFa} \\
            -f ${dstFa} \\
            -c ${chainFile} \\
            --reject ${unmVcf}

    # ---------- mini report ----------
    total=\$(bcftools view -H ${vcf} | wc -l)
    failed=\$(grep -vc '^#' ${unmVcf} || true)
    {
      echo "Total entries: \$total"
      echo "Failed to map: \$failed"
    } > ${id}_report.txt
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
