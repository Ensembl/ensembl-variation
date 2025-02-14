// extractMetadata parses and grabs the URN ID's entry from the pre-downloaded maveDB main.json file,
// which contains metadata for all URN ID's.
// It uses extract_meta.py to do this. 
process extract_metadata {
  // publishDir "${params.output}", mode: 'copy', overwrite: true

  input:
    val urn
    path metadata_file

  output:
    tuple val(urn), file("metadata.json")

  script:
  """
  echo "Processing URN: '${urn}'"
  python3 ${workflow.projectDir}/bin/extract_metadata.py --metadata_file ${metadata_file} --urn "${urn}"
  """
}