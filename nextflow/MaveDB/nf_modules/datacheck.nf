process datacheck_urns {
  tag "datacheck_urns"
  publishDir "reports", mode: 'copy', overwrite: true

  input:
    tuple path(final_gz), path(urn_file), path(log_csv)

  output:
    path("found_urns.txt")
    path("missing_urns.txt")
    path("missing_logs.csv")

  script:
  """
  #!/usr/bin/env bash
  set -euo pipefail

  # Extract URNs from the final output (gzip) and sort unique
  gzip -dc "${final_gz}" | awk -F'\\t' '
    NR==1 {
      for (i=1; i<=NF; i++) if (\$i=="urn") { col=i; break }
      if (!col) { print "ERROR: urn column not found in output header" > "/dev/stderr"; exit 2 }
      next
    }
    col && \$col && \$col!="urn" { print \$col }
  ' | sort -u > found_urns.txt

  # Sort unique input URNs (ignore blank lines)
  awk "NF>0" "${urn_file}" | sort -u > input_urns.txt

  # Compute missing URNs (present in input list but absent from output)
  comm -23 input_urns.txt found_urns.txt > missing_urns.txt

  # Collate logs for missing URNs (if any)
  if [[ -s missing_urns.txt && -s "${log_csv}" ]]; then
    (head -n1 "${log_csv}"; grep -F -f missing_urns.txt "${log_csv}") > missing_logs.csv || true
  else
    # Fallback header only
    echo '"file","time","urn","step","reason","subid"' > missing_logs.csv
  fi
  """
}
