#!/bin/bash

# # TEST
# urn="urn:mavedb:00000043-a-2"
# mappings_path="/nfs/production/flicek/ensembl/variation/jma/maveDB-test/mavedb_dbdump_data/mappings"
# scores_path="/nfs/production/flicek/ensembl/variation/jma/maveDB-test/mavedb_dbdump_data/scores"

echo "import_from_files - processing URN: '${urn}'" 2>&1

# Locate the mapping file using the original urn (with colons)
mapping_file=$(find ${mappings_path} -type f -iname "*${urn}*.json" | head -n 1)

# Replace colons with hyphens for searching the scores directory
score_urn=$(echo "${urn}" | sed 's/:/-/g')
score_file=$(find ${scores_path} -type f -iname "*${score_urn}.scores.csv" | head -n 1)

# Check if the mapping and scores files exist in the user-provided directories
if [ -z "${mapping_file}" ]; then
  echo "ERROR: No mapping file found for ${urn}" 2>&1
  exit 1
fi
if [ -z "${score_file}" ]; then
  echo "ERROR: No scores file found for ${urn}" 2>&1
  exit 1
fi

# If the files exist, report this
echo "Found mapping file: ${mapping_file} for ${urn}" 2>&1
echo "Found scores file: ${score_file} for ${urn} (using sanitized urn: ${score_urn})" 2>&1

# Copy mappings file to the working dir
cp "${mapping_file}" mappings.json

# Check if the file contains any lines starting with "tmp:"
if grep -q '^tmp:' ${score_file}; then

  echo "Score file for ${urn} contains temporary IDs (tmp:*). Replacing with file base name." 2>&1

  # Get the file's base name (e.g., "urn-mavedb-00000001-a-1")
  prefix=$(basename "${score_file}")
  prefix=${prefix%.scores.csv}
  echo "Using prefix: ${prefix}" 2>&1

  # Process the file with awk:
  # 1. Substitute "urn-mavedb-" with "urn:mavedb:"
  # 2. Print the columns
  # 3. For each subsequent row, substitute the first field: replace the pattern starting with "tmp:" up to and including the "#" with file prefix (file name) and "#"
  awk -F, -v OFS=, -v prefix="$prefix" '
    BEGIN {
      gsub(/^urn-mavedb-/, "urn:mavedb:", prefix)
    }
    NR==1 { print }
    NR>1 { sub(/^tmp:[^#]+#/, prefix "#" , $1); print }
  ' "${score_file}" >"scores.csv"

else
  # If the file does not contain any lines starting with "tmp:", copy the file as is
  echo "IDs as expected. Copying file as is." 2>&1
  cp "$score_file" ./scores.csv
fi

# If contents of pwd is mappings.json and scores.csv, then the files are copied successfully - check
if [ -f mappings.json ] && [ -f scores.csv ]; then
  echo "Files copied successfully" 2>&1
  echo "Contents of pwd: $(ls -l)" 2>&1
else
  echo "ERROR: Files not copied successfully" 2>&1
  echo "Contents of pwd: $(ls -l)" 2>&1
  exit 1
fi
