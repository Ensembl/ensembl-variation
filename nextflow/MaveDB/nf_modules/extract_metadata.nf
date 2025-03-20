// extractMetadata parses and grabs the URN ID's entry from the pre-downloaded maveDB main.json file,
// which contains metadata for all URN ID's.
// It then reformats it to match the JSON format expected by the rest of the pipeline which 
// was originally written to handle a different JSON structure. It uses extract_metadata.py to do this. 
process extract_metadata {
  // publishDir "${params.output}", mode: 'copy', overwrite: true

  input:
    val urn
    path metadata_file

  output:
    tuple val(urn), file("metadata.json"), file("LICENCE.txt")

  script:
  """
  #!/usr/bin/env bash
  
  set +e
  
  python3 ${workflow.projectDir}/bin/extract_metadata.py --metadata_file ${metadata_file} --urn "${urn}"
  
  # If metadata.json is missing or empty, create fallback files.
  if [ ! -s metadata.json ]; then
      echo "{}" > metadata.json
      echo "" > LICENCE.txt
      echo "WARNING: No metadata extracted for ${urn}, using fallback empty file." >&2
  fi
  """
}