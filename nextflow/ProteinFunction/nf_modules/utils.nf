#!/usr/bin/env nextflow

process decompress {
  /*
  Decompress file (if gzipped)

  Returns decompressed file
  */
  tag "${file}"
  container "quay.io/biocontainers/agat:0.9.0--pl5321hdfd78af_0"

  input:
    path file

  output:
    path '*', includeInputs: true

  """
  if [[ ${file.extension} == *gz ]]; then
    gunzip $file
  fi
  """
}
