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
include { get_hgvsp;
          run_variant_recoder;
          parse_vr_output } from './nf_modules/variant_recoder.nf'

include { map_scores_to_HGVSp_variants;
          map_scores_to_HGVSg_variants } from './nf_modules/mapping.nf'

log.info """
  Create MaveDB plugin data for VEP
  ---------------------------------
  mappings : ${params.mappings} 
  """

def split_by_mapping_type (files) {
  // split mapping files based on HGVS type (HGVSp or HGVSg files)

  tmp = files.map {
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
  mapping_files = Channel.fromPath( params.mappings + "/*.json" )
  mapping_files = split_by_mapping_type( mapping_files )
  
  // prepare HGVSp mappings
  get_hgvsp( mapping_files.hgvs_pro )
  run_variant_recoder( get_hgvsp.out )
  parse_vr_output( run_variant_recoder.out )
  map_scores_to_HGVSp_variants( parse_vr_output.out )

  // use MaveDB-prepared HGVSg mappings
  map_scores_to_HGVSg_variants( mapping_files.hgvs_nt )

  // lift-over coordinates to hg38 if needed

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
