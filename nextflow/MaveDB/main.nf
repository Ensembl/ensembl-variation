#!/usr/bin/env nextflow

/* 
 * Nextflow pipeline to create MaveDB data for VEP
 */

nextflow.enable.dsl=2

// Default params
params.help     = false
params.mappings = null

// Print usage
if (params.help) {
  log.info """
  Create MaveDB data for VEP
  --------------------------

  Usage:
    nextflow run main.nf -profile lsf -resume
  """
  exit 1
}

// Module imports
//include { decompress } from './nf_modules/utils.nf'

log.info """
  Crate MaveDB data for VEP
  -----------------------------
  mappings : ${params.mappings} 
  """

process split_by_mapping_type {
  tag "$file"

  input:
    path file

  output:
    tuple path(file), path('hgvsp.txt'), emit: hgvsp, optional: true
    tuple path(file), path('nucleotide.txt'), emit: input, optional: true

  """
  grep -m1 ":p." $file || cp $file nucleotide.txt

  if [[ ! -f nucleotide.txt ]]; then
    grep -Eo '".*:p..*"' $file | sed 's/"//g' | sed 's/value: //g' > hgvsp.txt
  fi

  # remove file if empty
  [[ -s hgvsp.txt ]] || rm -f hgvsp.txt
  """
}

process run_variant_recoder {
  tag "$file"

  input:
    tuple path(file), path(hgvs)

  output:
    tuple path(file), path('vr.json')

  """
  variant_recoder -i $hgvs --vcf_string > vr.json
  """
}

process prepare_vr_output {
  tag "$file"

  input:
    tuple path(file), path(vr)

  output:
    tuple path(file), path('vr.txt')

  """
  #!/usr/bin/env python3
  import csv
  import json

  f = open('$vr')
  data = json.load(f)

  line = []
  for result in data:
    for allele in result:
      info = result[allele]

      if type(info) is list and "Unable to parse" in info[0]:
        continue

      hgvsp  = info["input"]
      for string in info["vcf_string"]:
        chr, pos, ref, alt = string.split('-')
        line.append( [hgvsp, chr, pos, ref, alt] )
  f.close()

  w = open('vr.txt', 'w')
  writer = csv.writer(w, delimiter='\t')
  writer.writerows(line)
  w.close()
  """
}

workflow {
  mapping_files = Channel.fromPath(params.mappings + "/*")
  split_by_mapping_type(mapping_files)
  
  // prepare HGVSp mappings
  run_variant_recoder( split_by_mapping_type.out.hgvsp )
  prepare_vr_output( run_variant_recoder.out )

  // use HGVSg mappings
  // mix into single file
}

// Print summary
workflow.onComplete {
  println ( workflow.success ? """
        Workflow summary
        ----------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """ : """
        Failed: ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """
  )
}
