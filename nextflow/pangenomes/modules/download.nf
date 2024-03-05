process fetch_gene_symbol_lookup {
  // Download HGNC gene symbol lookup table with respective Ensembl identifiers

  output:
    path 'gene_symbol_table.txt', emit: file

  """
  # Download HGNC gene symbol table
  wget https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt
  awk -F"\t" '{if (\$20) print \$2"\t"\$20}' hgnc_complete_set.txt | awk 'NR > 1' > gene_symbol_table.txt
  """
}

process list_assemblies {
  // List available assemblies from URL (requires FTP protocol)
  input:
    val url
  output:
    stdout

  script:
    def link = url + "/"
  """
  curl -l ${link}
  """
}

process download_pangenomes_data {
  // Download pangenomes data

  input:
    val url
    val assembly
  output:
    tuple path('*.gtf.gz'), path('*.fa.gz')

  script:
    def link = url + "/" + assembly + "/"
  """
  wget -A "*genes.gtf.gz" --no-parent -r -nd ${link}
  wget -A "*unmasked.fa.gz" --no-parent -r -nd ${link}

  # remove older GTF files (based on alphabetical listing)
  ls *.gtf.gz | head -n -1 | xargs -r rm --
  """
}

