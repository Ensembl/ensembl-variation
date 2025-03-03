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

// Parameters for loading MaveDB from files:
params.from_files    = false       // Use local files instead of downloading via the MaveDB API
params.metadata_file = ""          // only used if from_files is true
params.mappings_path = ""          // only used if from_files is true
params.scores_path   = ""          // only used if from_files is true

// Print usage
if (params.help) {
  log.info """
  Create MaveDB data for VEP
  --------------------------

  Usage:
    ./main.nf -profile lsf -resume \
        --mappings [MaveDB_mappings_folder]

  Mandatory arguments:
    --urn           File with MaveDB URNs, such as 'urn:mavedb:00000001-a-1'

  Optional arguments:
    --ensembl       Path to Ensembl root directory (default: ${ENSEMBL_ROOT_DIR})
    --output        Path to output file (default: output/MaveDB_variants.tsv.gz)
    --registry      Path to Ensembl registry
    --from_files    Use local files instead of downloading via the MaveDB API
    --mappings_path Path to MaveDB mappings files (one JSON file per URN)
    --scores_path   Path to MaveDB scores files (one CSV file per URN)
    --metadata_file Path to MaveDB metadata file (one collated file, i.e. main.json)
    --licences      Comma-separated list of accepted licences (default: 'CC0')
    --round         Decimal places to round floats in MaveDB data (default: 4)
  """
  exit 1
}

// Module imports
include { filter_by_licence } from './subworkflows/filter.nf'
include { download_MaveDB_data } from './nf_modules/fetch.nf'
include { split_by_mapping_type } from './subworkflows/split.nf'
include { run_variant_recoder } from './nf_modules/variant_recoder.nf'
include { get_hgvsp } from './nf_modules/utils.nf'
include { map_scores_to_HGVSp_variants; map_scores_to_HGVSg_variants } from './nf_modules/mapping.nf'
include { download_chain_files; liftover_to_hg38 } from './nf_modules/liftover.nf'
include { concatenate_files; tabix } from './nf_modules/output.nf'
include { check_JVM_mem; print_params; print_summary } from '../utils/utils.nf'
include { import_from_files } from './nf_modules/import_from_files.nf'
include { extract_metadata } from './nf_modules/extract_metadata.nf'

// Main workflow
print_params('Create MaveDB plugin data for VEP', nullable=['registry'])
check_JVM_mem(min=50.4)
print_summary()

workflow {
  urn = Channel
      .fromPath(params.urn, checkIfExists: true)
      .splitText()
      .map { it.trim() }
      // .take(15)  // TEST: take the first n lines from the file

  // If --from_files is true, use local files instead of downloading via the MaveDB API
  if (params.from_files) {
    
    // metaChannel extracts the metadata from the large metadata file 
    // metaChannel outputs tuples of [urn, metadata.json, LICENSE.txt]
    metaChannel = extract_metadata(urn, params.metadata_file)

    // Log removed URNs (those that do NOT have a "CC0" license)
    metaChannel
        .filter { it[2].text != 'CC0' }
        .subscribe { println "NOTE: Discarded ${it[0]} based on license (${it[2].text})" }

    // Filter URNs based on license (only keep "CC0")
    filteredMetaChannel = metaChannel.filter { params.licences.contains(it[2].text) }

    // Remove LICENSE.txt, leaving output of filteredMetaChannel to be tuple: [urn, metadata.json]
    filteredMetaChannel = filteredMetaChannel.map { [it[0], it[1]] }

    // mainChannel finds and pulls in the mappings and scores files for each urn
    // mainChannel outputs tuples of [urn, mappings.json, scores.csv] - only with "CC0" license URNs
    mainChannel = import_from_files(filteredMetaChannel.map { it[0] })

    // Joins by URN ID and outputs tuple: [urn, mappings.json, scores.csv, metadata.json]
    files = mainChannel.join(filteredMetaChannel, by: 0).map { a, b, c, d -> return [
          urn      : a,  // URN ID
          mappings : b,  // Mappings file path
          scores   : c,  // Scores file path
          metadata : d   // Metadata file path
        ]}

  } else {
    // If --from_files is false, download MaveDB data via the API
    licences = params.licences.tokenize(",")
    urn = filter_by_licence(urn, licences)
    download_MaveDB_data(urn)
    files = download_MaveDB_data.out.map { [urn: it[0], mappings: it[1], scores: it[2], metadata: it[3]] }
  }
  
  // Split mapping.json files by mapping type - HGVSg or HGVSp
  files = split_by_mapping_type(files)

  // use MaveDB-prepared HGVSg mappings
  map_scores_to_HGVSg_variants(files.hgvs_nt)
  download_chain_files()
  liftover_to_hg38(map_scores_to_HGVSg_variants.out, download_chain_files.out)

  // prepare HGVSp mappings
  get_hgvsp(files.hgvs_pro)
  hgvsp = get_hgvsp.out.filter { it.last().size() > 0 }
  run_variant_recoder(hgvsp)
  map_scores_to_HGVSp_variants(run_variant_recoder.out)

  // concatenate output files into a single file
  output_files = liftover_to_hg38.out
                  .mix(map_scores_to_HGVSp_variants.out)
                  .collect { it.last() }
  concatenate_files(output_files)
  tabix(concatenate_files.out)
}
