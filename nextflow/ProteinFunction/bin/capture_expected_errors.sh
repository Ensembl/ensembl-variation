#!/bin/bash

# Example: capture_expected_errors.sh sift_align .command.err 54
md5=${1?:1st arg must be md5 of peptide sequence}
category=${2?:2nd arg must be the name of a set of errors}
stderr_file=${3?:3rd arg must be the error file}
exit_status=${4?:4th arg must be the exit status}
error_out=${5?:5th arg must be a file where to write errors}

# Set expected errors by category
if [ "${category}" == "sift_align" ]; then
  errors=("PSI-BLAST found no hits"
          "Not enough sequences found by the PSI-BLAST search"
          "Not enough sequences (only [0-9]*) found by the PSI-BLAST search")
elif [ "${category}" == "pph2" ]; then
  errors=("Failed to locate sequence position")
fi

# Capture expected errors
for error in "${errors[@]}"; do
  if grep -q "${error}" ${stderr_file}; then
    echo "======= Captured expected error: job is successful ======="
    cat ${stderr_file}
    echo "=========================================================="
    rm ${stderr_file}
    exit_status=0

    echo "${md5}\t${error}" > ${error_out}
    break
  fi
done

# If error is not captured, exit with error status
if [ ${exit_status} -ne 0 ]; then exit ${exit_status}; fi
