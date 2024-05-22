process download_chain_files {
  // Download UCSC LiftOver chain files

  output: path("*.over.chain.gz")

  shell:
  '''
  reference="Hg38"
  for genome in hg16 hg17 hg18 hg19; do
    url="https://hgdownload.cse.ucsc.edu/goldenpath/${genome}/liftOver/"
    name="${genome}To${reference}.over.chain.gz"
    
    wget ${url}/${name}
  done
  '''
}

process liftover_to_hg38 {
  // Identify genome of reference used to map variants and lift-over to hg38
  container 'quay.io/biocontainers/pyliftover:0.4--py_0'
  tag "${mappings.simpleName}"

  input: 
    tuple path(mappings), path(mapped_variants)
    path(chain_files)
  output:
    tuple path(mappings), path("liftover_*.tsv")

  """
  liftover.py --mappings ${mappings} \
              --mapped_variants ${mapped_variants} \
              --reference hg38
  """
}
