#!/usr/bin/env nextflow

/* 
 * Nextflow pipeline to remap VCF files
 */

nextflow.enable.dsl=2
description = "Pipeline to remap VCF files"
separator   = "-" * description.length()

// Default params
params.help   = false
params.ids    = null
params.vcf    = null

params.chain_dir = null
params.fasta_dir = null

params.out_dir = 'output'
params.report = 'crossmap_report.txt'

// Print usage
if (params.help) {
  log.info """
  ${description}
  ${separator}

  Usage:
    nextflow run -profile slurm -resume main.nf \\
      --ids ids.txt \\
      --vcf sample.vcf.gz \\
      --chain_dir /path/to/chain/dir \\
      --fasta_dir /path/to/fasta/dir

  Mandatory arguments:
    --vcf        VCF file to remap
    --ids        File containing IDs of genomes
    --chain_dir  Path to chain directory (filenames must have the ID given in --ids)
    --fasta_dir  Path to fasta directory (filenames must have the ID given in --ids)

  Optional arguments:
    --out_dir    Path to output directory (default: 'output')
    --report     Filename with remapping statistics (default: 'crossmap_report.txt')
  """
  exit 1
}

include { crossmap; tabix; report } from './modules/crossmap.nf'
include { print_params; print_summary } from '../utils/utils.nf'

print_params(description, separator)
print_summary()

workflow {
  ids = Channel.fromPath(params.ids, checkIfExists: true).splitText { it.trim() }

  vcf       = Channel.fromPath(params.vcf, checkIfExists: true)
  chain_dir = Channel.fromPath(params.chain_dir, checkIfExists: true)
  fasta_dir = Channel.fromPath(params.fasta_dir, checkIfExists: true)

  crossmap(vcf, ids, chain_dir, fasta_dir)
  tabix(crossmap.out.vcf)
  report(crossmap.out.report.collect())
}
