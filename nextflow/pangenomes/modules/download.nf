process fetch_gene_symbol_lookup {
  // Download HGNC gene symbol lookup table with respective Ensembl identifiers

  output:
    path 'gene_symbol_table.txt', emit: file

  """
  # Download HGNC gene symbol table
  wget https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt
  awk -F"\t" '{if (\$20) print \$2"\t"\$20}' hgnc_complete_set.txt | awk 'NR > 1' > gene_symbol_table.txt
  """
}

process list_assemblies {
  // List available assemblies from URL
  input:
    val url
  output:
    path 'assemblies.txt'

  script:
    def link = url.replaceAll('/+$', '') + '/'
  """
  curl --fail --silent --show-error --location ${link} \
    | sed -n 's#.*href="\\(GCA_[0-9][0-9]*\\.[0-9][0-9]*\\)/".*#\\1#p' \
    | sort -u > assemblies.txt
  test -s assemblies.txt
  """
}

process download_pangenomes_data {
  // Download pangenomes data
  tag "${assembly}"

  input:
    val url
    val assembly
  output:
    tuple val(assembly), path('genes.gff3.gz'), path('unmasked.fa.gz')

  script:
    def link = url.replaceAll('/+$', '') + '/' + assembly
  """
  release=\$(curl --fail --silent --show-error --location ${link}/ensembl/geneset/ \
    | sed -n 's#.*href="\\([0-9][0-9][0-9][0-9]_[0-9][0-9]\\)/".*#\\1#p' \
    | sort -V | tail -n 1)
  test -n "\${release}"
  wget -O genes.gff3.gz ${link}/ensembl/geneset/\${release}/genes.gff3.gz
  wget -O unmasked.fa.gz ${link}/genome/unmasked.fa.gz
  gzip -t genes.gff3.gz
  gzip -t unmasked.fa.gz
  """
}
