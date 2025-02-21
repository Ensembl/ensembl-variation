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
params.from_files    = false                                                // default: use local files
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
    --from_files    Use local files instead of downloading from MaveDB

  Optional arguments:
    --ensembl       Path to Ensembl root directory (default: ${ENSEMBL_ROOT_DIR})
    --output        Path to output file (default: output/MaveDB_variants.tsv.gz)
    --registry      Path to Ensembl registry
    --mappings_path Path to MaveDB mappings files (per URN)
    --scores_path   Path to MaveDB scores files (per URN)
    --metadata_file Path to MaveDB metadata files (main.json, one collated file)

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
include { map_scores_to_HGVSp_variants;
      map_scores_to_HGVSg_variants } from './nf_modules/mapping.nf'
include { download_chain_files;
      liftover_to_hg38 } from './nf_modules/liftover.nf'
include { concatenate_files;
      tabix } from './nf_modules/output.nf'
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
      .take(2)  // TEST: take the first n lines from the file

  // Choose which module to use based on --from_files (true/false)
  if (params.from_files) {
    // metaChannel outputs tuples of [urn, metadata.json]
    metaChannel = extract_metadata(urn, params.metadata_file)
    
    // mainChannel outputs tuples of [urn, mappings.json, scores.csv]
    mainChannel = import_from_files(urn)

    // Joins by URN ID and outputs tuples of [urn, mappings.json, scores.csv, metadata.json]
    files = mainChannel.join(metaChannel, by: 0).map { a, b, c, d -> return [
          urn      : a,  // URN ID
          mappings : b,  // Mappings file path
          scores   : c,  // Scores file path
          metadata : d   // Metadata file path
        ]}

  } else {
    // Download_MaveDB_data using the API
    licences = params.licences.tokenize(",")
    urn = filter_by_licence(urn, licences)
    download_MaveDB_data(urn)
    files = download_MaveDB_data.out.map { [urn: it[0], mappings: it[1], scores: it[2], metadata: it[3]] }
  }

  // // Split mapping files based on HGVS type (HGVSp or HGVSg files) - make 2 channels 
  // files_split = files.map { entry -> 
  //         def hgvs = "unknown"
  //         def mappings_file = file(entry.mappings)

  //         mappings_file.withReader { reader ->
  //             reader.eachLine { line ->
  //                 if (line.contains("hgvs.")) {
  //                     hgvs = line.contains("hgvs.p") ? "hgvs.p" : "hgvs.g"
  //                     return
  //                 }
  //             }
  //         }
          
  //         return entry + [hgvs: hgvs]
  //     }
  //     .branch {
  //         hgvs_pro: it.hgvs == "hgvs.p",
  //         hgvs_nt:  it.hgvs == "hgvs.g"
  //     }

  files = split_by_mapping_type(files)

  println("_pro")
  files.hgvs_pro.view()
  println("_nt")
  files.hgvs_nt.view()

  // // use MaveDB-prepared HGVSg mappings
  map_scores_to_HGVSg_variants(files.hgvs_nt)
  download_chain_files()
  liftover_to_hg38(map_scores_to_HGVSg_variants.out, download_chain_files.out)

  // // prepare HGVSp mappings
  get_hgvsp(files.hgvs_pro)
  hgvsp = get_hgvsp.out.filter { it.last().size() > 0 }
  run_variant_recoder(hgvsp)
  map_scores_to_HGVSp_variants(run_variant_recoder.out)

  // // concatenate output files into a single file
  // output_files = liftover_to_hg38.out
  //                 .mix(map_scores_to_HGVSp_variants.out)
  //                 .collect { it.last() }
  // concatenate_files(output_files)
  // tabix(concatenate_files.out)
}
