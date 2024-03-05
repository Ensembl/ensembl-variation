process create_latest_annotation {
  // Use VEP to create Phenotypes or GO annotation for latest assembly

  input:
    val plugin
    val version
    val species
    val user
    val host
    val port
  output:
    path "*.gz", emit: file
    path "*.gz.tbi", emit: tbi

  afterScript 'rm *.tmp'

  script:
  def plugin_opts = plugin == 'GO' ?
    'GO,match=gene_symbol' :
    'Phenotypes,dir=.,include_types=Gene'
  """
  ${params.vep} --database --db_version ${version} \
      --species $species \
      --user ${user} --host ${host} --port ${port} \
      --plugin ${plugin_opts} \
    || echo "Avoid VEP error: Cannot detect format from STDIN"
  """
}

process filter_Phenotypes_gene_annotation {
  // Filter Phenotypes annotation to only contain Gene data

  input:
    path annotation
  output:
    path annotation.baseName, emit: file

  """
  zcat ${annotation} | awk -F" " '\$3 == "Gene"' > ${annotation.baseName}
  """
}

process create_pangenomes_annotation {
  container 'docker://biocontainers/pandas:1.5.1_cv1'

  input:
    val plugin
    val version
    tuple path(gtf), path(fasta)
    path annotation
    path gene_symbols

  output:
    tuple path(gtf), path(fasta), path('*.g*f')

  script:
    def opts = (plugin == 'GO' ? '--go' : '--pheno') + " ${annotation}"
    def lookup = (gene_symbols.name == 'null') ? '' : "--gene_symbols ${gene_symbols}"
  """
  create_pangenomes_annotation.py --version ${version} --gtf ${gtf} ${opts} ${lookup}
  """
}

process tabix_plugin_annotation {
  publishDir "${params.outdir}", mode: 'copy', pattern: '*plugin*.gz*'

  input:
    tuple path(gtf), path(fasta), path(annotation)
  output:
    tuple path(gtf), path(fasta), path('*.gz'), path('*.gz.tbi')

  """
  bgzip ${annotation}
  tabix ${annotation}.gz
  """
}
