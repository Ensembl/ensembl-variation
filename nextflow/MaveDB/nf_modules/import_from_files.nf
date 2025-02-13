// File: nf_modules/import_from_files.nf

process import_from_files {
  // This process now only receives a simple urn value.
  input:
    val urn

  // It outputs a tuple: (urn, mappings file, scores file) - need to add metadata file later, once extraction of URN ID is implemented.
  output:
    tuple val(urn), path("mappings.json"), path("scores.csv")
    
  script:
  """
  # Locate the mapping file using the original urn (with colons)
  mapping_file=\$(find ${params.mappings_path} -type f -iname "*${urn}*.json" | head -n 1)
  
  # Replace colons with hyphens for searching the scores directory.
  score_urn=\$(echo ${urn} | sed 's/:/-/g')
  score_file=\$(find ${params.scores_path} -type f -iname "*\${score_urn}*.csv" -not -iname "*counts*" | head -n 1)
  
  if [ -z "\$mapping_file" ]; then
      echo "ERROR: No mapping file found for \$urn" >&2
      exit 1
  fi
  if [ -z "\$score_file" ]; then
      echo "ERROR: No scores file found for \$urn" >&2
      exit 1
  fi
  
  echo "Found mapping file: \$mapping_file for \$urn"
  echo "Found scores file: \$score_file for \$urn (using sanitized urn: \$score_urn)"
  
  cp "\$mapping_file" mappings.json
  cp "\$score_file" scores.csv
  
  echo "Files successfully copied for \$urn"
  """
}