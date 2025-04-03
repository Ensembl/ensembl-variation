process map_scores_to_HGVSp_variants {
  // Download MaveDB scores and map associated variants by HGVSp

  tag "${urn}"
  input:  tuple val(urn), path(mappings), path(scores), path(metadata), path(vr)
  output: tuple val(urn), path('map_*.tsv')

  memory { mappings.size() * 4.B + 1.GB }

  script:
  def round = params.round ? "--round ${params.round}" : ""

  // If --from_files is true, use local files instead of downloading via the MaveDB API
  def script_name = params.from_files ? "map_scores_to_variants_fromfiles.py" : "map_scores_to_variants.py"

  """
  #!/usr/bin/env bash
  
  set +e
  
  ${script_name} --urn ${urn} \\
                 --scores ${scores} \\
                 --mappings ${mappings} \\
                 --metadata ${metadata} \\
                 --vr $vr \\
                 ${round} \\
                 --output map_${urn}.tsv

  # Check if the output file exists and is non-empty, if not, create an empty file
  if [ ! -s map_${urn}.tsv ]; then
      echo "WARNING: map_${urn}.tsv is empty or doesn't exist. Creating fallback empty file." >&2
      echo "" > map_${urn}.tsv
  fi
  """
}

process map_scores_to_HGVSg_variants {
  // Download MaveDB scores and map associated variants from HGVSg

  tag "${urn}"
  input:  tuple val(urn), path(mappings), path(scores), path(metadata), val(hgvs)
  output: tuple val(urn), path(metadata), path('*map_*.tsv')

  memory { mappings.size() * 2.B + 1.GB }

  script:
  def round = params.round ? "--round ${params.round}" : ""

  // If --from_files is true, use local files instead of downloading via the MaveDB API
  def script_name = params.from_files ? "map_scores_to_variants_fromfiles.py" : "map_scores_to_variants.py"

  """
  set +e

  ${script_name} --urn ${urn} \\
                 --scores ${scores} \\
                 --mappings ${mappings} \\
                 --metadata ${metadata} \\
                 ${round} \\
                 --output map_${urn}.tsv

  # Check if the output file exists and is non-empty, if not, create an empty file
  if [ ! -s map_${urn}.tsv ]; then
      echo "WARNING: map_${urn}.tsv is empty or doesn't exist. Creating fallback empty file." >&2
      echo "" > map_${urn}.tsv
  fi
  """
}