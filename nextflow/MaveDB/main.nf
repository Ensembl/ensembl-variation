#!/usr/bin/env nextflow

/* 
 * Nextflow pipeline to create MaveDB plugin data for VEP
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
include { run_variant_recoder; prepare_vr_mappings } from './nf_modules/variant_recoder.nf'

log.info """
  Create MaveDB plugin data for VEP
  ---------------------------------
  mappings : ${params.mappings} 
  """

def split_by_mapping_type (files) {
  // split mapping files based on HGVS type (HGVSp or HGVSg files)

  tmp = Channel.fromPath(files).map {
    it.withReader {
      while( line = it.readLine() ) {
        if (line.contains("hgvs.")) {
          // get first line describing HGVS type
          if (line.contains("hgvs.p")) {
            hgvs = "hgvs.p"
          } else if (line.contains("hgvs.g")) {
            hgvs = "hgvs.g"
          } else {
            throw new Exception("Error: HGVS type in '${line.trim()}' not expected")
          }
          break
        }
      }
    }
    [file: it, hgvs: hgvs]
  }.branch{
    hgvs_pro: it.hgvs == "hgvs.p"
    hgvs_nt:  it.hgvs == "hgvs.g"
  }
  
  // clean up: only return the files in each branch
  files = [hgvs_pro: null, hgvs_nt: null]
  files.hgvs_pro = tmp.hgvs_pro.map { it.file }
  files.hgvs_nt  = tmp.hgvs_nt.map { it.file }
  return files
}

workflow {
  mapping_files = Channel.fromPath( params.mappings + "/*" )
  mapping_files = split_by_mapping_type( mapping_files )
  
  // prepare HGVSp mappings
  run_variant_recoder( mapping_files.hgvs_pro )
  prepare_vr_mappings( run_variant_recoder.out )

  // use MaveDB-prepared HGVSg mappings
  process_MaveDB_mappings( mapping_files.hgvs_nt )

  // concatenate data into a single file
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
