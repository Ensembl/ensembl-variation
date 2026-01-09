process collate_logs {
  tag "collate_logs"
  publishDir "reports", mode: 'copy', overwrite: true

  input:
    tuple val(workdir), val(trigger)

  output:
    path("collated_logs.csv")

  script:
  """
  #!/usr/bin/env bash
  set -euo pipefail
  
  collate_logs.sh "${workdir}" collated_logs.csv
  """
}
