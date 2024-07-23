#!/usr/bin/env nextflow

/* 
 * Nextflow pipeline to create MaveDB plugin data for VEP
 */

nextflow.enable.dsl=2

// Default params
params.help     = false
params.urn      = null
params.ensembl  = "${ENSEMBL_ROOT_DIR}"
params.output   = "output/MaveDB_variants.tsv.gz"
params.registry = null

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
    --urn       File with MaveDB URNs, such as 'urn:mavedb:00000001-a-1'

  Optional arguments:
    --ensembl   Path to Ensembl root directory (default: ${ENSEMBL_ROOT_DIR})
    --output    Path to output file (default: output/MaveDB_variants.tsv.gz)
    --registry  Path to Ensembl registry

    --licences  Comma-separated list of accepted licences (default: 'CC0')
    --round     Decimal places to round floats in MaveDB data (default: 4)
  """
  exit 1
}

// Module imports
include { filter_by_licence } from './subworkflows/filter.nf'
include { download_MaveDB_data } from './nf_modules/fetch.nf'
include { split_by_mapping_type } from './subworkflows/split.nf'
include { run_variant_recoder } from './nf_modules/variant_recoder.nf'
include { get_hgvsp } from './nf_modules/utils.nf'
include { map_scores_to_HGVSp_variants;
          map_scores_to_HGVSg_variants } from './nf_modules/mapping.nf'
include { download_chain_files;
          liftover_to_hg38 } from './nf_modules/liftover.nf'
include { concatenate_files;
          tabix } from './nf_modules/output.nf'

log.info """
  Create MaveDB plugin data for VEP
  ---------------------------------
  urn      : ${params.urn}
  output   : ${params.output}
  ensembl  : ${params.ensembl}
  registry : ${params.registry}

  licences : ${params.licences}
  round    : ${params.round}
  """

workflow {
  // filter MaveDB URNs based on file-specific licence
  urn = Channel
          .fromPath( params.urn, checkIfExists: true )
          .splitText()
          .map { it.trim() }
  licences = params.licences.tokenize(",")
  urn = filter_by_licence(urn, licences)

  // download mappings and scores from MaveDB
  download_MaveDB_data(urn)
  files = download_MaveDB_data.out
            .map { [urn: it[0], mappings: it[1], scores: it[2], metadata: it[3]] }
  files = split_by_mapping_type( files )

  // use MaveDB-prepared HGVSg mappings
  map_scores_to_HGVSg_variants( files.hgvs_nt )
  download_chain_files()
  liftover_to_hg38( map_scores_to_HGVSg_variants.out,
                    download_chain_files.out )

  // prepare HGVSp mappings
  get_hgvsp( files.hgvs_pro )
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
