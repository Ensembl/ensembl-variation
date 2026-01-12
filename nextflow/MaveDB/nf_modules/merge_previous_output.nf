process merge_previous_output {
  tag "merge_previous_output"

  input:
    tuple path(current_combined), path(prev_output), path(skipped_urns)

  output:
    path("merged_combined.tsv")

  script:
  """
  #!/usr/bin/env bash
  set -euo pipefail

  if [[ ! -f "${current_combined}" ]]; then
    echo "ERROR: current combined file not found: ${current_combined}" >&2
    exit 1
  fi

  if [[ ! -f "${prev_output}" ]]; then
    echo "ERROR: previous output file not found: ${prev_output}" >&2
    exit 1
  fi

  urn_col=\$(head -n1 "${current_combined}" | awk -F'\\t' '{for(i=1;i<=NF;i++) if(\$i=="urn"){print i; exit}}')
  if [[ -z "\${urn_col}" ]]; then
    echo "ERROR: urn column not found in ${current_combined}" >&2
    exit 2
  fi

  # Decompress previous output once
  gzip -dc "${prev_output}" > prev.tsv

  # Extract only rows belonging to skipped URNs (if any)
  if [[ -s "${skipped_urns}" ]]; then
    awk -F'\\t' -v OFS='\\t' -v col="\${urn_col}" 'NR==FNR {u[\$0]=1; next} NR==1 {print; next} (\$col in u)' "${skipped_urns}" prev.tsv > prev_subset.tsv
  else
    : > prev_subset.tsv
  fi

  # Merge: header once, then current rows, then recovered previous rows for skipped URNs
  {
    head -n1 "${current_combined}"
    tail -n +2 "${current_combined}"
    if [[ -s prev_subset.tsv ]]; then
      tail -n +2 prev_subset.tsv
    fi
  } | awk 'NR==1 {print; next} !seen[\$0]++' > merged_combined.tsv
  """
}
