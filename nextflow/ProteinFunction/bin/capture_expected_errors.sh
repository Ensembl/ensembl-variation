#!/bin/bash

# Example: capture_expected_errors.sh 37e493a9feb018c2ce69e68c0093f6d2 sift_align .command.error 54 expected_error.txt
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
elif [ "${category}" == "sift" ]; then
  errors=("The following sequences have been removed because they  were found to be over 100% identical with your protein query")
elif [ "${category}" == "pph2" ]; then
  errors=("Failed to locate sequence position")
fi

# catch memory error
mem_error="Some of the step tasks have been OOM Killed."
if grep -q "${mem_error}" ${stderr_file}; then
  exit 140
fi

# Capture expected errors
for error in "${errors[@]}"; do
  if grep -q "${error}" ${stderr_file}; then
    echo "======= Captured expected error: job is successful ======="
    cat ${stderr_file}
    echo "=========================================================="
    rm ${stderr_file}
    exit_status=0

    echo "${md5}\t${error}\t${category}" > ${error_out}
    break
  fi
done

# If error is not captured, exit with error status
if [ ${exit_status} -ne 0 ]; then exit ${exit_status}; fi
