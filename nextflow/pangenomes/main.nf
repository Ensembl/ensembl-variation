#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.url     = 'ftp://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/'
params.outdir  = 'outdir'
params.sw      = "${ENSEMBL_ROOT_DIR}"
params.vep     = "${params.sw}/ensembl-vep/vep"

params.version = null
params.species = "homo_sapiens"

params.user    = null
params.host    = null
params.port    = null

log.info "\nCreate GO and Phenotype annotations for pangenomes"
log.info "=================================================="
for (a in params) {
  // print param
  log.info "  ${a.getKey().padRight(8)} : ${a.getValue()}"
  // raise error if param is null
  if (!a.getValue()) exit 1, "ERROR: parameter --${a.getKey()} not defined"  
}
log.info ""

include { fetch_gene_symbol_lookup;
          list_assemblies; 
          download_pangenomes_data } from './modules/download.nf'
include { create_latest_annotation;
          filter_Phenotypes_gene_annotation;
          tabix_plugin_annotation;
          create_pangenomes_annotation } from './modules/annotation.nf'
include { tabix_gtf; decompress_fasta } from './modules/utils.nf'
include { test_annotation } from './modules/test.nf'

workflow create_go_annotations {
  take:
    data
  main:
    go_grch38 = create_latest_annotation('GO', params.version, params.species,
                                         params.user, params.host, params.port)
    go_pan = create_pangenomes_annotation('GO', params.version, data, go_grch38.file, '/')
    go_pan = go_pan | tabix_plugin_annotation | decompress_fasta | tabix_gtf
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
    pheno_pan = pheno_pan | tabix_plugin_annotation | decompress_fasta | tabix_gtf
    test_annotation('Phenotypes', pheno_pan)
}

workflow {
  assemblies = list_assemblies(params.url).splitCsv().flatten()
  pangenomes = download_pangenomes_data(params.url, assemblies).transpose()
  create_go_annotations(pangenomes)
  create_pheno_annotations(pangenomes)
}
