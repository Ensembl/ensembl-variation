
process import_from_files {
  // Receives an urn value from the main workflow
  input:
    val urn

  // Outputs a tuple: (urn, mappings file, scores file) - need to add metadata file later, once extraction of URN ID is implemented.
  output:
    tuple val(urn), file("mappings.json"), file("scores.csv")
    
  script:
  """
  echo "import_from_files - processing URN: '\${urn}'" 2>&1
  
  # Locate the mapping file using the original urn (with colons)
  mapping_file=\$(find \${params.mappings_path} -type f -iname "*\${urn}*.json" | head -n 1)
  
  # Replace colons with hyphens for searching the scores directory
  score_urn=\$(echo "\${urn}" | sed 's/:/-/g')
  score_file=\$(find \${params.scores_path} -type f -iname "*\${score_urn}*scores.csv" | head -n 1)
  
  if [ -z "\$mapping_file" ]; then
      echo "ERROR: No mapping file found for \${urn}" 2>&1
      exit 1
  fi
  if [ -z "\$score_file" ]; then
      echo "ERROR: No scores file found for \${urn}" 2>&1
      exit 1
  fi
  
  echo "Found mapping file: \$mapping_file for \${urn}" 2>&1
  echo "Found scores file: \$score_file for \${urn} (using sanitized urn: \$score_urn)" 2>&1
  
  # Copy mappings file to the working dir with standardised name
  cp "\$mapping_file" mappings.json

  # Check if the file contains any lines starting with "tmp:"
  if grep -q '^tmp:' \$score_file; then

    echo "Score file for \${urn} contains tmp: IDs. Replacing with file base name"  2>&1

    # Get the file's base name (e.g., "urn-mavedb-00000001-a-1")
    prefix=$(basename \$score_file .scores.csv)

    # Process the file with awk:
    # 1. Substitute "urn-mavedb-" with "urn:mavedb:"
    # 2. Print the columns 
    # 3. For each subsequent row, substitute the first field: replace the pattern starting with "tmp:" up to and including the "#" with file prefix (file name) and "#"
    awk -F, -v OFS=, -v prefix="\$prefix" '
      BEGIN { 
        gsub(/^urn-mavedb-/, "urn:mavedb:", prefix)
      }
      NR==1 { print }
      NR>1 { sub(/^tmp:[^#]+#/, prefix"#", $1); print }
    ' \$score_file > "\$prefix.scores.csv"

  fi
  
  echo "Mapping and Score files successfully copied for ${urn}"  2>&1
  """
}