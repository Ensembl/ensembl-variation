process map_scores_to_HGVSp_variants {
  // Download MaveDB scores and map associated variants by HGVSp

  tag "${mappings.simpleName}"
  input:  tuple path(mappings), path(vr)
  output: tuple path(mappings), path('map_*.tsv')

  script:
  def urn = "${mappings.simpleName}"
  """
  map_scores_to_variants.py --urn ${urn} \
                            --mappings ${mappings} \
                            --matches $vr \
                            --output map_${urn}.tsv
  """
}

process map_scores_to_HGVSg_variants {
  // Download MaveDB scores and map associated variants from HGVSg

  tag "${mappings.simpleName}"
  input:  path(mappings)
  output: tuple path(mappings), path('map_*.tsv')

  script:
  def urn = "${mappings.simpleName}"
  """
  map_scores_to_variants.py --urn ${urn} \
                            --mappings ${mappings} \
                            --output map_${urn}.tsv
  """
}
