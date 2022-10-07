/* 
 * Nextflow pipeline to predict protein function using SIFT and PolyPhen-2
 */

nextflow.enable.dsl=2

// Default params
params.help       = false
params.outdir     = "outdir"
params.gtf        = null
params.fasta      = null
params.translated = null
params.species    = "homo_sapiens"

// Variation database params
params.host     = null
params.port     = null
params.user     = null
params.pass     = null
params.database = null

// SIFT params
params.sift_run_type = "NONE"
params.median_cutoff = 2.75 // as indicated in SIFT's README
params.blastdb       = null

// PolyPhen-2 params
params.pph_run_type = "NONE"
params.pph_data     = null

// print usage
if (params.help) {
  log.info """
  Predict protein function using SIFT and PolyPhen-2
  --------------------------------------------------

  Usage:
    nextflow run main.nf -profile lsf -resume \
             --gtf basenji.gtf,boxer.gtf --fasta basenji.fa,boxer.fa \
             --pph_run_type  UPDATE --pph_data [path/to/pph_data] \
             --sift_run_type UPDATE --blastdb  [path/to/blastdb] \
             --host [h] --port [p] --user [u] --pass [p] --database [db]

  General options:
    --gtf FILE             Comma-separated annotation GTF files; requires FASTA
    --fasta FILE           Comma-separated FASTA files with genomic sequences;
                           requires GTF
    --translated FILE      Comma-separated FASTA files with peptide sequence;
                           skips sequence translation based on GTF and FASTA
    --outdir DIRNAME       Name of output dir (default: "outdir")
    --species VAL          Latin species name (default: homo_sapiens);
                           PolyPhen-2 only works for human

  Database options (mandatory):
    --host VALUE           Database server host
    --port VALUE           Database server port
    --user VALUE           Database server user
    --pass VALUE           Database server password
    --database VALUE       Database name

  SIFT options:
    --sift_run_type VALUE  SIFT run type:
                             - "FULL" to run for all translations
                             - "UPDATE" to run for new/changed translations
                             - "NONE" to exclude this analysis (default)
    --blastdb DIR          SIFT-formatted BLAST database directory
                           (e.g., uniref100)
    --median_cutoff VALUE  Protein alignment's median cutoff (default: 2.75)

  PolyPhen-2 options:
    --pph_run_type VALUE   PolyPhen-2 run type:
                             - "FULL" to run for all translations
                             - "UPDATE" to run for new/changed translations
                             - "NONE" to exclude this analysis (default)
    --pph_data DIR         Path to PolyPhen-2 data
  """
  exit 1
}

log.info """
  Predict protein function using SIFT and PolyPhen-2
  --------------------------------------------------
  GTF       : ${params.gtf}
  FASTA     : ${params.fasta}

  host      : ${params.host}
  port      : ${params.port}
  user      : ${params.user}
  database  : ${params.database}

  SIFT run  : ${params.sift_run_type}
  blastdb   : ${params.blastdb}

  PPH2 run  : ${params.pph_run_type}
  PPH2 data : ${params.pph_data}
  """

// Module imports
include { translate_fasta }           from './nf_modules/translations.nf'
include { store_translation_mapping } from './nf_modules/database_utils.nf'
include { run_sift_pipeline }         from './nf_modules/sift.nf'
include { run_pph2_pipeline }         from './nf_modules/polyphen2.nf'

// Check run type for each protein function predictor
def check_run_type ( run ) {
  run_types = ["NONE", "UPDATE", "FULL"]
  if ( !run_types.contains( run ) ) {
    supported = run_types.join(', ')
    exit 1, "ERROR: invalid run type $run! Supported run types are $supported"
  }
}
check_run_type( params.sift_run_type )
check_run_type( params.pph_run_type )

// Check if supplying PolyPhen-2 data and if species is human 
if ( params.pph_run_type != "NONE" ) {
  if ( params.species != "homo_sapiens" ) {
    exit 1, "ERROR: PolyPhen-2 only works with human data"
  } else if ( !params.pph_data) {
    exit 1, "ERROR: --pph_data must be supplied when running PolyPhen-2"
  }
}

// Check blastdb for SIFT
if ( params.blastdb ) {
  params.blastdb_name = file(params.blastdb).name
  params.blastdb_dir = file(params.blastdb).parent
} else if ( params.sift_run_type != "NONE" ) {
  exit 1, "ERROR: --blastdb must be supplied when running SIFT"
}

workflow {
  // Translate transcripts from GTF and FASTA if no translation FASTA is given
  if (!params.translated) {
    gtf   = Channel.fromPath(  params.gtf.tokenize(','), checkIfExists: true)
    fasta = Channel.fromPath(params.fasta.tokenize(','), checkIfExists: true)
    translate_fasta(gtf, fasta)
    translated = translate_fasta.out
  } else {
    translated = Channel.fromPath(params.translated.tokenize(','))
  }

  // Parse translation FASTA file
  translated = translated
                 .splitFasta(record: [id: true, seqString: true, text: true])
                 .unique() // remove duplicated entries
                 .map{ it -> [id: it.id,
                              // Remove stop codon (asterisk)
                              text: it.text.replaceAll(/\*/, ""),
                              seqString: it.seqString.replaceAll(/\*/, ""),
                              // Add MD5 hashes for the translation sequences
                              md5: it.seqString.replaceAll(/\*/, "").md5() ]}

  // Write translation mapping with transcript ID and MD5 hashes to database
  translation_mapping = translated.collectFile(
                          name: "translation_mapping.tsv",
                          storeDir: params.outdir,
                          newLine: true) { it.id + "\t" + it.md5 }
  store_translation_mapping(translation_mapping)

  // Get unique translations based on MD5 hashes of their sequences
  translated = translated.unique { it.md5 }

  // Run protein function prediction
  if ( params.sift_run_type != "NONE" ) run_sift_pipeline( translated )
  if ( params.pph_run_type  != "NONE" ) run_pph2_pipeline( translated )
}

// Print summary
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
