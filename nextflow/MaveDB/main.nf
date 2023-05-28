#!/usr/bin/env nextflow

/* 
 * Nextflow pipeline to create MaveDB plugin data for VEP
 */

nextflow.enable.dsl=2

// Default params
params.help     = false
params.mappings = null
params.round    = 4

// Print usage
if (params.help) {
  log.info """
  Create MaveDB data for VEP
  --------------------------

  Usage:
    ./main.nf -profile lsf -resume \
              --mappings [MaveDB_mappings_folder]

  Mandatory arguments:
    --mappings  Path to directory containing MaveDB mapping files in JSON format

  Optional arguments:
    --round     Decimal places to round floats in MaveDB data (default: 4)

  Note: the mappings must be named with the MaveDB accession, such as
        '00000001-a-1.json'.
  """
  exit 1
}

// Module imports
include { get_hgvsp;
          run_variant_recoder } from './nf_modules/variant_recoder.nf'
include { check_if_open_access;
          map_scores_to_HGVSp_variants;
          map_scores_to_HGVSg_variants } from './nf_modules/mapping.nf'
include { download_chain_files;
          liftover_to_hg38 } from './nf_modules/liftover.nf'

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

process round_scores {
  // Round MaveDB scores to given number of decimal places

  input:
    path file
    val precision
  output:
    path "rounded_scores.txt"

  """
  col=\$(awk -v RS='\t' '/score/{print NR; exit}' ${file})
  awk 'NR>1 {\$12=sprintf("%.${precision}f", \$12)}1' ${file} > rounded_scores.txt
  """
}

process tabix {
  input:  path out
  output: path "${name}"

  script:
  def name = "MaveDB_variants.tsv.gz"
  """
  sort -k1 -nk2,3 ${out} > ${name.baseName}
  bgzip ${name.baseName}
  tabix -s1 -b2 -e3 ${name}
  """
}

workflow {
  check_if_open_access( Channel.fromPath( params.mappings + "/*.json" ) )
  mapping_files = split_by_mapping_type( check_if_open_access.out )

  // use MaveDB-prepared HGVSg mappings
  map_scores_to_HGVSg_variants( mapping_files.hgvs_nt )
  download_chain_files()
  liftover_to_hg38( map_scores_to_HGVSg_variants.out,
                    download_chain_files.out )

  // prepare HGVSp mappings
  get_hgvsp( mapping_files.hgvs_pro )
  run_variant_recoder( get_hgvsp.out )
  map_scores_to_HGVSp_variants( run_variant_recoder.out )

  // concatenate output files into a single file
  output_files = liftover_to_hg38.out.
                   mix(map_scores_to_HGVSp_variants.out).
                   collect { it.last() }
  concatenate_files( output_files )
  tabix( concatenate_files.out )
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
