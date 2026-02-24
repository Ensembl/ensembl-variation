#!/usr/bin/env nextflow

import java.nio.file.Files
import java.nio.file.StandardCopyOption

/*
 * Nextflow pipeline to create MaveDB plugin data for VEP
 */

nextflow.enable.dsl=2

// Default params
params.help     = false
params.urn      = null
params.previous_urn = null
params.previous_output = null
params.ensembl  = "${ENSEMBL_ROOT_DIR}"
params.output   = "output/MaveDB_variants.tsv.gz"
params.registry = null

params.licences = "CC0" // Open-access only
params.round    = 4

// Parameters for loading MaveDB from files:
params.from_files    = true        // Use local files instead of downloading via the MaveDB API
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
    --from_files    Use local files instead of downloading via the MaveDB API (default: true, this is advised)
    --previous_urn  Path to a previous URN list; any URN already present there is skipped
    --previous_output Path to a previous run's output (.tsv.gz); skipped URNs are taken from here and merged back
    --mappings_path Path to MaveDB mappings files (one JSON file per URN)
    --scores_path   Path to MaveDB scores files (one CSV file per URN)
    --metadata_file Path to MaveDB metadata file (one collated file, i.e. main.json)
    --licences      Comma-separated list of accepted licences (default: 'CC0')
    --round         Decimal places to round floats in MaveDB data (default: 4)
  """
  exit 1
}

def read_urn_file(path) {
  def urn_file = new File(path.toString())
  if (!urn_file.exists()) {
    exit 1, "ERROR: URN file not found: ${urn_file}"
  }

  return urn_file.readLines()
                    .collect { it.trim() }
                    .findAll { it }
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
include { collate_logs } from './nf_modules/collate_logs.nf'
include { datacheck_urns } from './nf_modules/datacheck.nf'
include { merge_previous_output } from './nf_modules/merge_previous_output.nf'

// Main workflow
print_params('Create MaveDB plugin data for VEP', nullable=['registry', 'previous_urn', 'previous_output'])
check_JVM_mem(min=50.4)
print_summary()

workflow {
  def urn_list = read_urn_file(params.urn)
  if (params.previous_output && !file(params.previous_output).exists()) {
    exit 1, "ERROR: previous output not found: ${params.previous_output}"
  }
  def previous_urns = params.previous_urn ? read_urn_file(params.previous_urn).toSet() : [] as Set
  def filtered_urns = previous_urns ? urn_list.findAll { !previous_urns.contains(it) } : urn_list
  def skipped_urns = urn_list - filtered_urns

  if (params.previous_urn) {
    def skipped = urn_list.size() - filtered_urns.size()
    log.info "Loaded ${previous_urns.size()} URNs from ${params.previous_urn}; skipping ${skipped} already processed URNs; ${filtered_urns.size()} remaining."
  }

  if (!filtered_urns) {
    if (params.previous_output) {
      log.info "No new URNs to process; reusing ${params.previous_output} as final output"
      Files.copy(file(params.previous_output).toPath(), file(params.output).toPath(), StandardCopyOption.REPLACE_EXISTING)
      def prevIndex = file(params.previous_output + ".tbi")
      if (prevIndex.exists()) {
        Files.copy(prevIndex.toPath(), file(params.output + ".tbi").toPath(), StandardCopyOption.REPLACE_EXISTING)
      }
      exit 0
    } else {
      exit 0, "No URNs to process after applying --previous_urn filter"
    }
  }

  def skipped_urns_file = workflow.workDir.resolve('skipped_urns.txt')
  skipped_urns_file.toFile().parentFile.mkdirs()
  skipped_urns_file.toFile().text = skipped_urns ? skipped_urns.join(System.lineSeparator()) + System.lineSeparator() : ""

  def urn_file_for_datacheck
  if (params.previous_output) {
    urn_file_for_datacheck = file(params.urn) // expect final merged output to contain all URNs
  } else if (params.previous_urn) {
    urn_file_for_datacheck = workflow.workDir.resolve('urns_to_process.txt')
    urn_file_for_datacheck.toFile().parentFile.mkdirs()
    urn_file_for_datacheck.toFile().text = filtered_urns.join(System.lineSeparator()) + System.lineSeparator()
  } else {
    urn_file_for_datacheck = file(params.urn)
  }

  urn = Channel.from(filtered_urns)

  // If --from_files is true, use local files instead of downloading via the MaveDB API
  if (params.from_files) {
    
    // metaChannel extracts the metadata from the large metadata file 
    // metaChannel outputs tuples of [urn, metadata.json, LICENSE.txt]
    metaChannel = extract_metadata(urn, params.metadata_file)

    // Log removed URNs (those that do NOT have a "CC0" license)
    metaChannel
    // Log removed URNs (e.g. - those that do NOT have a "CC0" license)
    metaChannel
        .filter { !params.licences.tokenize(",").contains(it[2].text) }
        .subscribe { println "NOTE: Discarded ${it[0]} based on license (${it[2].text})" }

    // Filter URNs based on license
    filteredMetaChannel = metaChannel.filter { params.licences.tokenize(",").contains(it[2].text) }

    // Remove LICENSE.txt, leaving output of filteredMetaChannel to be tuple: [urn, metadata.json]
    filteredMetaChannel = filteredMetaChannel.map { [it[0], it[1]] }

    // mainChannel finds and pulls in the mappings and scores files for each urn
    // mainChannel outputs tuples of [urn, mappings.json, scores.csv]
    mainChannel = import_from_files(filteredMetaChannel.map { it[0] })

    // Joins by URN ID and outputs tuple: [urn, mappings.json, scores.csv, metadata.json]
    files = mainChannel.join(filteredMetaChannel, by: 0).map { a, b, c, d -> return [
          urn      : a,  // URN ID
          mappings : b,  // Mappings JSON file path
          scores   : c,  // Scores CSV file path
          metadata : d   // Metadata JSON file path
        ]}

  } else {
    // If --from_files is false, download MaveDB data via the API (this option is not advised as it's unreliable)
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

  def tabix_input = concatenate_files.out
  if (params.previous_output) {
    tabix_input = merge_previous_output(tabix_input.map { combined -> tuple(combined, file(params.previous_output), file(skipped_urns_file)) })
  }

  def tabix_out = tabix(tabix_input)

  // collate logs into a single csv
  log_collation_input = tabix_out.map { t -> [workflow.workDir.toString(), t] }
  collate_logs(log_collation_input)

  // datacheck: compare final URNs against input list and pull logs for any missing URNs
  datacheck_input = collate_logs.out.map { logs -> tuple(file(params.output), urn_file_for_datacheck, logs) }
  datacheck_urns(datacheck_input)
}
