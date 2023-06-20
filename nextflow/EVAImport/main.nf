/* 
 * Nextflow pipeline to import variants from EVA
 */

nextflow.enable.dsl=2
nextflow.enable.strict = true

// Scripts
eva_script         = "${ENSEMBL_ROOT_DIR}/ensembl-variation/scripts/import/import_vcf.pl"
var_syn_script     = "${ENSEMBL_ROOT_DIR}/ensembl-variation/scripts/import/import_variant_synonyms"
var_set_script     = "${ENSEMBL_ROOT_DIR}/ensembl-variation/scripts/import/import_set_from_file.pl"
var_set_script_2   = "${ENSEMBL_ROOT_DIR}/ensembl-variation/scripts/import/post_process_variation_feature_variation_set.pl"
copy_tables_script = "${ENSEMBL_ROOT_DIR}/ensembl-internal-variation/scripts/copy_tables_after_eva.sh"
citations_script   = "${ENSEMBL_ROOT_DIR}/ensembl-internal-variation/scripts/import_citation_EVA.pl"

// Common params
params.help            = false
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
params.old_host        = ""
params.old_port        = ""
params.old_dbname      = ""

// Params to copy tables after import
params.host            = ""
params.dbname          = ""

// Params for citations
params.citations_file  = null

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

// Print usage
if (params.help) {
  log.info """
  ------------------------------
  Import variation data from EVA
  ------------------------------
  
  Usage:
    nextflow run main.nf \\
             --species sus_scrofa \\
             --registry ensembl.registry \\
             --release 111 \\
             --input_file GCA_000003025.6_current_ids.vcf.gz \\
             --var_syn_file GCA_000003025.6_merged_ids.vcf.gz \\
             --output_file output_import.log \\
             --host [new database host] \\
             --dbname [new database name] \\
             --old_dbname [previous database name] \\

    Options (mandatory):
    --species             species name
    --registry            registry file pointing to variation and core databases
    --release             release number
    --input_file          EVA input file
    --var_syn_file        EVA variation synonyms file [GCA_*_merged_ids.vcf.gz]
    --output_file         output file to write number of skipped variants in the EVA import
    --host                new variation database host (necessary to prepare the db for the import)
    --dbname              new variation database name (necessary to prepare the db for the import)
    --old_dbname          previous variation database name (necessary to prepare the db for the import)

    Options (only mandatory for rat):
    --old_host             previous variation database host
    --old_port             previous variation database port

    Options (optional):
    --chr_synonyms        file that contains the chromosome synonyms
    --citations_file      text file with list of rs ids linked to publication id from the previous database
  
  """
 exit 1
}

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

if(!params.old_dbname) {
  exit 1, "ERROR: please provide the previous database name to copy phenotype and SV tables from the previous database"
}

if( (!params.old_host || !params.old_port || !params.old_dbname) && params.species == "rattus_norvegicus") {
  exit 1, "ERROR: please provide a host (--old_host), port (--old_port) and db name (--old_dbname) for a previous rat database"
}

if(!params.output_file) {
  exit 1, "ERROR: please provide an output file to the EVA import script"
}

// Build command to run EVA import
command_to_run = " -i ${params.input_file} --source ${params.source} --source_description '${params.description}' --version ${params.version} --registry ${params.registry} --species ${params.species} --skip_tables '${params.skip_tables}'"

log.info """
  Importing EVA data with the following parameters: \
  ${command_to_run}
"""

process run_eva {
  cpus "${params.fork}"

  input:
  val wait
  path eva_script
  val options
  val merge_all_types
  val fork
  val sort_vf
  val chr_synonyms
  val remove_prefix
  val output_file
  
  output: val 'ok'
  
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
  val wait
  path var_syn_script
  val source_name
  val species
  val input_file
  val registry
  val host
  val port
  val dbname
  
  output: val 'ok'
  
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
  val wait
  val my_species_set
  val var_set_script
  val files_path
  val filenames
  val species
  val registry

  output: val 'ok'

  script:

  def input_file = filenames.get(my_species_set)
  """
  perl ${var_set_script} -load_file ${files_path}${input_file} -registry ${registry} -species ${species} -variation_set ${my_species_set}
  """
}

// Post-process variation feature sets
// Populates variation_feature.variation_set_id
process run_variation_set_2 {
  input:
  val wait
  val var_set_script_2
  val species
  val registry

  output: val 'ok'

  script:

  """
  perl ${var_set_script_2} -registry_file ${registry} -species ${species}
  """
}

process prepare_tables {
  input:
  val copy_script
  val host
  val dbname
  val old_dbname

  output: val 'ok'

  script:
  def new_host = "${host}-ensadmin"

  """
  sh ${copy_script} -h ${new_host} -o ${old_dbname} -n ${dbname}
  """
}

process run_citations {
  input:
  val wait
  path citations_script
  val species
  val registry
  val file

  output: val 'ok'

  script:

  """
  perl ${citations_script} -load_file ${file} -registry ${registry} -species ${species}
  """
}


workflow {
  prepare_tables(copy_tables_script, params.host, params.dbname, params.old_dbname)

  run_eva(prepare_tables.out, file(eva_script), command_to_run, params.merge_all_types, params.fork, params.sort_vf, params.chr_synonyms, params.remove_prefix, params.output_file)

  run_variant_synonyms(run_eva.out, file(var_syn_script), params.source, params.species, params.var_syn_file, params.registry, params.old_host, params.old_port, params.old_dbname)

  if(set_names[params.species]) {
    my_species_set = Channel.fromList(set_names.get(params.species))
    run_variation_set(run_variant_synonyms.out, my_species_set, var_set_script, files_path, filenames, params.species, params.registry)

    run_variation_set_2(run_variation_set.out.collect(), var_set_script_2, params.species, params.registry)
  }

  if(params.citations_file) {
    run_citations(run_variant_synonyms.out, file(citations_script), params.species, params.registry, params.citations_file)
  }
}

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