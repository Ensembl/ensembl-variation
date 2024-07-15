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
params.keep_id   = false

params.out_dir = 'output'
params.report = 'crossmap_report.txt'

lookup = [
  "HG02257.1"    : "GCA_018466835.1",
  "HG02257.2"    : "GCA_018466845.1",
  "HG02559.1"    : "GCA_018466855.1",
  "HG02559.2"    : "GCA_018466985.1",
  "HG02486.1"    : "GCA_018467005.1",
  "HG02486.2"    : "GCA_018467015.1",
  "HG01891.2"    : "GCA_018467155.1",
  "HG01891.1"    : "GCA_018467165.1",
  "HG01258.2"    : "GCA_018469405.1",
  "HG03516.1"    : "GCA_018469415.1",
  "HG03516.2"    : "GCA_018469425.1",
  "HG01123.2"    : "GCA_018469665.1",
  "HG01258.1"    : "GCA_018469675.1",
  "HG01361.2"    : "GCA_018469685.1",
  "HG01123.1"    : "GCA_018469695.1",
  "HG01361.1"    : "GCA_018469705.1",
  "HG01358.2"    : "GCA_018469865.1",
  "HG02622.2"    : "GCA_018469875.1",
  "HG02622.1"    : "GCA_018469925.1",
  "HG02717.2"    : "GCA_018469935.1",
  "HG02630.1"    : "GCA_018469945.1",
  "HG02630.2"    : "GCA_018469955.1",
  "HG01358.1"    : "GCA_018469965.1",
  "HG02717.1"    : "GCA_018470425.1",
  "HG02572.1"    : "GCA_018470435.1",
  "HG02572.2"    : "GCA_018470445.1",
  "HG02886.2"    : "GCA_018470455.1",
  "HG02886.1"    : "GCA_018470465.1",
  "HG01175.1"    : "GCA_018471065.1",
  "HG01106.1"    : "GCA_018471075.1",
  "HG01175.2"    : "GCA_018471085.1",
  "HG00741.2"    : "GCA_018471095.1",
  "HG00741.1"    : "GCA_018471105.1",
  "HG01106.2"    : "GCA_018471345.1",
  "HG00438.2"    : "GCA_018471515.1",
  "HG02148.1"    : "GCA_018471525.1",
  "HG02148.2"    : "GCA_018471535.1",
  "HG01952.2"    : "GCA_018471545.1",
  "HG01952.1"    : "GCA_018471555.1",
  "HG00673.2"    : "GCA_018472565.1",
  "HG00621.1"    : "GCA_018472575.1",
  "HG00673.1"    : "GCA_018472585.1",
  "HG00438.1"    : "GCA_018472595.1",
  "HG00621.2"    : "GCA_018472605.1",
  "HG01071.2"    : "GCA_018472685.1",
  "HG01928.2"    : "GCA_018472695.1",
  "HG01928.1"    : "GCA_018472705.1",
  "HG00735.1"    : "GCA_018472715.1",
  "HG01071.1"    : "GCA_018472725.1",
  "HG00735.2"    : "GCA_018472765.1",
  "HG03579.2"    : "GCA_018472825.1",
  "HG03579.1"    : "GCA_018472835.1",
  "HG01978.1"    : "GCA_018472845.1",
  "HG03453.2"    : "GCA_018472855.1",
  "HG01978.2"    : "GCA_018472865.1",
  "HG03540.2"    : "GCA_018473295.1",
  "HG03453.1"    : "GCA_018473305.1",
  "HG03540.1"    : "GCA_018473315.1",
  "HG03486.1"    : "GCA_018503245.1",
  "NA18906.2"    : "GCA_018503255.1",
  "NA18906.1"    : "GCA_018503285.1",
  "HG03486.2"    : "GCA_018503525.1",
  "HG02818.1"    : "GCA_018503575.1",
  "HG02818.2"    : "GCA_018503585.1",
  "HG01243.1"    : "GCA_018504045.1",
  "HG02080.1"    : "GCA_018504055.1",
  "HG02723.2"    : "GCA_018504065.1",
  "HG02723.1"    : "GCA_018504075.1",
  "HG02080.2"    : "GCA_018504085.1",
  "HG01109.2"    : "GCA_018504365.1",
  "HG01243.2"    : "GCA_018504375.1",
  "NA20129.1"    : "GCA_018504625.1",
  "NA20129.2"    : "GCA_018504635.1",
  "HG01109.1"    : "GCA_018504645.1",
  "NA21309.2"    : "GCA_018504655.1",
  "NA21309.1"    : "GCA_018504665.1",
  "HG02109.2"    : "GCA_018505825.1",
  "HG03492.1"    : "GCA_018505835.1",
  "HG03492.2"    : "GCA_018505845.1",
  "HG02055.1"    : "GCA_018505855.1",
  "HG02109.1"    : "GCA_018505865.1",
  "HG02055.2"    : "GCA_018506125.1",
  "HG03098.1"    : "GCA_018506155.1",
  "HG03098.2"    : "GCA_018506165.1",
  "HG00733.1"    : "GCA_018506955.1",
  "HG00733.2"    : "GCA_018506975.1",
  "HG02145.2"    : "GCA_018852585.1",
  "HG02145.1"    : "GCA_018852595.1",
  "T2T_CHM13_v2" : "GCA_009914755.4"
]

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
    --keep_id    Boolean: rename filename based on hard-coded lookup table (default) or keep original ID?
  """
  exit 1
}

include { crossmap; rename_and_tabix; report } from './modules/crossmap.nf'
include { print_params; print_summary } from '../utils/utils.nf'

print_params(description, separator)
print_summary()

workflow {
  ids = Channel.fromPath(params.ids, checkIfExists: true).splitText { it.trim() }

  vcf       = Channel.fromPath(params.vcf, checkIfExists: true)
  chain_dir = Channel.fromPath(params.chain_dir, checkIfExists: true)
  fasta_dir = Channel.fromPath(params.fasta_dir, checkIfExists: true)

  crossmap(vcf, ids, chain_dir, fasta_dir)
  rename_and_tabix(crossmap.out.vcf, lookup)
  report(crossmap.out.report.collect(), lookup)
}
