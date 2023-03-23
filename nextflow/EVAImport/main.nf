/* 
 * Nextflow pipeline to import variants from EVA
 */

nextflow.enable.dsl=2
nextflow.enable.strict = true

// Scripts
eva_script     = "${ENSEMBL_ROOT_DIR}/ensembl-variation/scripts/import/import_vcf.pl"
var_syn_script = "${ENSEMBL_ROOT_DIR}/ensembl-variation/scripts/import/import_variant_synonyms"

// Common params
//params.help            = false
params.species         = null
params.release         = null
params.registry        = null

// Params for EVA import script
params.input_file      = null
params.source          = "EVA"
params.description     = "Short variant data imported from EVA"
params.version         = 4
params.remove_prefix   = false
params.chr_synonyms    = false
params.merge_all_types = true
params.fork            = 10
params.skip_tables     = "allele,allele_code,population,population_genotype,genotype_code,compressed_genotype_var,sample"
params.output_file     = null
params.sort_vf         = true

// Params for variant synonyms import
params.var_syn_file    = null


// Check input params
if(!params.species) {
  exit 1, "ERROR: species name must be provided when running EVA import"
}

if(!params.input_file || !file(params.input_file)) {
  exit 1, "ERROR: a valid input file must be provided when running EVA import"
}

if(!params.release || !params.registry) {
  exit 1, "ERROR: release version and registry file must be provided when running EVA import"
}


// Build command to run
command_to_run = " -i ${params.input_file} --source ${params.source} --source_description '${params.description}' --version ${params.version} --registry ${params.registry} --species ${params.species} --skip_tables '${params.skip_tables}'"

if(params.merge_all_types) {
command_to_run += " --merge_all_types"
}

if(params.fork) {
command_to_run += " --fork ${params.fork}"
}

if(params.chr_synonyms) {
command_to_run += " --chr_synonyms ${params.chr_synonyms}"
}

if(params.remove_prefix) {
command_to_run += " --remove_prefix ${params.remove_prefix}"
}

if(params.sort_vf) {
command_to_run += " --sort_vf"
}

log.info """
  Import EVA script: ${eva_script} \
  Options: ${command_to_run}
"""


process run_eva {
  input:
  path eva_script
  val options
  val output_file
  
  output:
  
  script:
  """
  perl ${eva_script} ${options} --output_file ${output_file}
  """
}

process run_variant_synonyms {
  input:
  path var_syn_script
  val source_name
  val species
  val input_file
  val registry
  
  output:
  
  script:
  
  if(species == "sus_scrofa")
      """
        perl ${var_syn_script} --source_name ${source_name} --species ${species} --data_file ${input_file} --registry ${registry}
        perl ${var_syn_script} --source_name "pig_chip" --species ${species} --registry ${registry}
      """
  
  else if(species == "rattus_norvegicus")
      """
        perl ${var_syn_script} --source_name ${source_name} --species ${species} --data_file ${input_file} --registry ${registry}
        perl ${var_syn_script} --source_name "rat" --species ${species} --registry ${registry}
      """
  else 
      """
        perl ${var_syn_script} --source_name ${source_name} --species ${species} --data_file ${input_file} --registry ${registry}
      """
}


workflow {
  // TODO: run script to truncate tables

  run_eva(file(eva_script), command_to_run, params.output_file)
  run_variant_synonyms(file(var_syn_script), params.source, params.species, params.var_syn_file, params.registry)
  
}
