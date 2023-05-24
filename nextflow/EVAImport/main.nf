/* 
 * Nextflow pipeline to import variants from EVA
 */

nextflow.enable.dsl=2
nextflow.enable.strict = true

// Scripts
eva_script     = "${ENSEMBL_ROOT_DIR}/ensembl-variation/scripts/import/import_vcf.pl"
var_syn_script = "${ENSEMBL_ROOT_DIR}/ensembl-variation/scripts/import/import_variant_synonyms"
var_set_script = "${ENSEMBL_ROOT_DIR}/ensembl-variation/scripts/import/import_set_from_file.pl"

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
params.host            = ""
params.port            = ""
params.dbname          = ""

// Params for sets import
files_path = "/nfs/production/flicek/ensembl/variation/data/genotyping_chips/"
filenames  = [ "MGP":"mouse/mgp_set/mgp_variation_set.txt.gz",
               "PorcineHD":"Pig/GGP_Porcine_HD_ids.txt.gz",
               "PorcineLD":"Pig/GGP_Porcine_LD_ids.txt.gz",
               "PorcineSNP60":"Pig/Illumina_PorcineSNP60_ids.txt.gz",
               "Affy_PorcineHD":"Pig/Axiom_PigHD_v1_ids.txt.gz",
               "OvineSNP50":"Sheep_Illumina/OvineSNP50_ids.txt.gz",
               "OvineHDSNP":"Sheep_Illumina/OvineHDSNP_ids.txt.gz",
               "Chicken600K":"Chicken/Chicken600K_ids.txt.gz",
               "Illumina_EquineSNP50":"Horse/EquineSNP50_ids.txt.gz",
               "GoatSNP50":"Goat/GoatSNP50_ids.txt.gz",
               "BovineHD":"Cow/BovineHD_ids.txt.gz",
               "BovineLD":"Cow/BovineLD_C_ids.txt",
               "BovineSNP50":"Cow/BovineSNP50_ids.txt.gz"
             ]
set_names  = [ "mus_musculus":["MGP"],
              "sus_scrofa":["PorcineHD", "PorcineLD", "PorcineSNP60", "Affy_PorcineHD"],
              "ovis_aries":["OvineSNP50", "OvineHDSNP"],
              "gallus_gallus":["Chicken600K"],
              "equus_caballus":["Illumina_EquineSNP50"],
              "capra_hircus":["GoatSNP50"],
              "bos_taurus":["BovineHD", "BovineLD", "BovineSNP50"]
             ]


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

if( (!params.host || !params.port || !params.dbname) && params.species == "rattus_norvegicus") {
  exit 1, "ERROR: please provide a host, port and db name for a previous rat database"
}


// Build command to run EVA import
command_to_run = " -i ${params.input_file} --source ${params.source} --source_description '${params.description}' --version ${params.version} --registry ${params.registry} --species ${params.species} --skip_tables '${params.skip_tables}'"


log.info """
  Import EVA script: ${eva_script} \
  Options: ${command_to_run}
"""


process run_eva {
  input:
  path eva_script
  val options
  val merge_all_types
  val fork
  val sort_vf
  val chr_synonyms
  val remove_prefix
  val output_file
  
  output:
  
  script:
  def sort_vf_table     = sort_vf ? " --sort_vf" : ""
  def merge_all         = merge_all_types ? " --merge_all_types" : ""
  def use_fork          = fork ? "--fork ${fork}" : ""
  def chr_synonyms_file = chr_synonyms ? " --chr_synonyms ${chr_synonyms}" : ""
  def rm_prefix         = remove_prefix ? " --remove_prefix ${remove_prefix}" : ""
  
  """
  perl ${eva_script} ${options} $sort_vf_table $merge_all $use_fork $chr_synonyms_file $rm_prefix --output_file ${output_file}
  """
}

process run_variant_synonyms {
  input:
  path var_syn_script
  val source_name
  val species
  val input_file
  val registry
  val host
  val port
  val dbname
  
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
        perl ${var_syn_script} --source_name "rat" --species ${species} --registry ${registry} --host ${host} --port ${port} --user 'ensro' --db_name $dbname
      """
  else 
      """
        perl ${var_syn_script} --source_name ${source_name} --species ${species} --data_file ${input_file} --registry ${registry}
      """
}

process run_variation_set {
  input:
  path var_set_script
  val files_path
  val filenames
  val set_names
  val species
  val registry

  exec:
  
  // variable my_set has to be before any if statement
  // related to this issue: https://github.com/nextflow-io/nextflow/issues/804
  def my_set = set_names.get(species)
  for (String name : my_set) {
    def input_file = filenames.get(name)
    """
      perl ${var_set_script} -load_file ${files_path}${input_file} -registry ${registry} -species ${species} -variation_set ${name}
    """
  }

}


workflow {
  // TODO: run script to truncate tables

  run_eva(file(eva_script), command_to_run, params.merge_all_types, params.fork, params.sort_vf, params.chr_synonyms, params.remove_prefix, params.output_file)
  run_variant_synonyms(file(var_syn_script), params.source, params.species, params.var_syn_file, params.registry, params.host, params.port, params.dbname)
  
  // variation_set has to be populated before import
  if(set_names[params.species]) {
    run_variation_set(file(var_set_script), files_path, filenames, set_names, params.species, params.registry)
  }
}
