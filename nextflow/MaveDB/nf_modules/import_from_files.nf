process import_from_files {
  // Receives an urn value from the main workflow
  input:
    val urn

  // Outputs a tuple: (urn, mappings file, scores file) - need to add metadata file later, once extraction of URN ID is implemented.
  output:
    tuple val(urn), file("mappings.json"), file("scores.csv")
    
  script:
  """
  set +e
  
  # Export environment variables so bash script can access them
  export urn=${urn}
  export mappings_path=${params.mappings_path}
  export scores_path=${params.scores_path}
  
  # Call bash import script
  bash ${workflow.projectDir}/bin/import_from_files.sh

  # Check if the mappings.json file exists and is non-empty, if not, create an empty file
  if [ ! -s mappings.json ]; then
      echo "{}" > mappings.json
      echo "WARNING: mappings.json is empty for ${urn}; using fallback empty JSON." >&2
  fi
  
  # Check if the scores.csv file exists and is non-empty, if not, create an empty file
  if [ ! -s scores.csv ]; then
      echo "" > scores.csv
      echo "WARNING: scores.csv is empty for ${urn}; using fallback empty file." >&2
  fi

  echo "Import from files completed for ${urn}"
  """
}