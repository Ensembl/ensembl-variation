#!/usr/bin/env nextflow

/* 
 * Nextflow pipeline to create MaveDB plugin data for VEP
 */

nextflow.enable.dsl=2

// Default params
params.help     = false
params.mappings = null
params.output   = "./MaveDB_variants.tsv.gz"
params.licences = "CC0" // Open-access only
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
    --output    Path to output file (default: './MaveDB_variants.tsv.gz')
    --licences  Comma-separated list of accepted licences (default: 'CC0')
    --round     Decimal places to round floats in MaveDB data (default: 4)

  Note: the mappings must be named with the MaveDB accession, such as
        '00000001-a-1.json'.
  """
  exit 1
}

// Module imports
include { get_files_by_licence;
          split_by_mapping_type } from './nf_modules/licences.nf'
include { get_hgvsp;
          run_variant_recoder } from './nf_modules/variant_recoder.nf'
include { map_scores_to_HGVSp_variants;
          map_scores_to_HGVSg_variants } from './nf_modules/mapping.nf'
include { download_chain_files;
          liftover_to_hg38 } from './nf_modules/liftover.nf'
include { concatenate_files;
          tabix } from './nf_modules/output.nf'

log.info """
  Create MaveDB plugin data for VEP
  ---------------------------------
  mappings : ${params.mappings}
  output   : ${params.output}

  licences : ${params.licences}
  round    : ${params.round}
  """

workflow {
  // prepare data based on licence
  files = Channel.fromPath( params.mappings + "/*.json" )
  licences = params.licences.tokenize(",")
  filtered = get_files_by_licence(files, licences)
  mapping_files = split_by_mapping_type( filtered )

  // use MaveDB-prepared HGVSg mappings
  map_scores_to_HGVSg_variants( mapping_files.hgvs_nt )
  download_chain_files()
  liftover_to_hg38( map_scores_to_HGVSg_variants.out,
                    download_chain_files.out )

  // prepare HGVSp mappings
  get_hgvsp( mapping_files.hgvs_pro )
  hgvsp = get_hgvsp.out.filter { it.last().size() > 0 }
  run_variant_recoder( hgvsp )
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
