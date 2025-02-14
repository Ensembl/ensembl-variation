
process import_from_files {
  // Receives an urn value from the main workflow
  input:
    val urn

  // Outputs a tuple: (urn, mappings file, scores file) - need to add metadata file later, once extraction of URN ID is implemented.
  output:
    tuple val(urn), file("mappings.json"), file("scores.csv")
    
  script:
  """
  echo "Processing URN: '${urn}'"
  
  # Locate the mapping file using the original urn (with colons)
  mapping_file=\$(find ${params.mappings_path} -type f -iname "*${urn}*.json" | head -n 1)
  
  # Replace colons with hyphens for searching the scores directory
  score_urn=\$(echo "${urn}" | sed 's/:/-/g')
  score_file=\$(find ${params.scores_path} -type f -iname "*\${score_urn}*scores.csv" | head -n 1)
  
  if [ -z "\$mapping_file" ]; then
      echo "ERROR: No mapping file found for ${urn}" >&2
      exit 1
  fi
  if [ -z "\$score_file" ]; then
      echo "ERROR: No scores file found for ${urn}" >&2
      exit 1
  fi
  
  echo "Found mapping file: \$mapping_file for ${urn}"
  echo "Found scores file: \$score_file for ${urn} (using sanitized urn: \$score_urn)"
  
  # Create symbolic links or copy files to the working directory.
  # Here we use cp to create local files named mappings.json and scores.csv.
  cp "\$mapping_file" mappings.json
  cp "\$score_file" scores.csv
  
  echo "Mapping and Score files successfully copied for ${urn}"
  """
}