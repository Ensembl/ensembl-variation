#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.url     = 'https://ftp.ebi.ac.uk/pub/ensemblorganisms/Homo_sapiens/'
params.input_manifest = null
params.outdir  = 'outdir'
params.sw      = "${ENSEMBL_ROOT_DIR}"
params.vep     = "${params.sw}/ensembl-vep/vep"

params.version = null
params.species = "homo_sapiens"

params.user    = null
params.host    = null
params.port    = null

include { fetch_gene_symbol_lookup;
          list_assemblies;
          download_pangenomes_data } from './modules/download.nf'
include { create_latest_annotation;
          filter_Phenotypes_gene_annotation;
          tabix_plugin_annotation;
          create_pangenomes_annotation } from './modules/annotation.nf'
include { tabix_gff3; decompress_fasta } from './modules/utils.nf'
include { test_annotation } from './modules/test.nf'

include { check_JVM_mem; print_params; print_summary } from '../utils/utils.nf'
print_params('Create GO and Phenotype annotations for pangenomes',
             nullable=['input_manifest'])
check_JVM_mem(min=0.4)
print_summary()

workflow create_go_annotations {
  take:
    data
  main:
    go_grch38 = create_latest_annotation('GO', params.version, params.species,
                                         params.user, params.host, params.port)
    go_pan = create_pangenomes_annotation('GO', params.version, data, go_grch38.file, '/')
    go_pan = go_pan | tabix_plugin_annotation | decompress_fasta | tabix_gff3
    test_annotation('GO', go_pan)
}

workflow create_pheno_annotations {
  take:
    data
  main:
    pheno_grch38 = create_latest_annotation('Phenotypes', params.version, params.species,
                                            params.user, params.host, params.port)
    filtered     = filter_Phenotypes_gene_annotation(pheno_grch38.file)
    gene_symbol  = fetch_gene_symbol_lookup()
    pheno_pan    = create_pangenomes_annotation('Phenotypes', params.version,
                                                data, filtered, gene_symbol.file)
    pheno_pan = pheno_pan | tabix_plugin_annotation | decompress_fasta | tabix_gff3
    test_annotation('Phenotypes', pheno_pan)
}

workflow {
  if (params.input_manifest) {
    pangenomes = Channel
      .fromPath(params.input_manifest, checkIfExists: true)
      .splitCsv(header: true, sep: '\t')
      .map { row ->
        tuple(row.accession,
              file(row.gff3, checkIfExists: true),
              file(row.fasta, checkIfExists: true))
      }
  }
  else {
    assemblies = list_assemblies(params.url)
      .splitText()
      .map { it.trim() }
      .filter { it }
    pangenomes = download_pangenomes_data(params.url, assemblies)
  }
  create_go_annotations(pangenomes)
  create_pheno_annotations(pangenomes)
}
