process import_from_files {
  // Receives an urn value from the main workflow
  input:
    val urn

  // Outputs a tuple: (urn, mappings file, scores file) - need to add metadata file later, once extraction of URN ID is implemented.
  output:
    tuple val(urn), file("mappings.json"), file("scores.csv")
    
  script:
  """
  # Export environment variables so bash script can access them
  export urn=${urn}
  export mappings_path=${params.mappings_path}
  export scores_path=${params.scores_path}
  
  # Call bash script
  bash ${workflow.projectDir}/bin/import_from_files.sh
  """
}