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

// Print usage
if (params.help) {
  log.info """
  Predict protein function using SIFT and PolyPhen-2
  --------------------------------------------------

  Usage:
    nextflow run main.nf -profile lsf -resume \\
             --species canis_lupus_familiaris \\
             --gtf basenji.gtf,boxer.gtf --fasta basenji.fa,boxer.fa \\
             --pph_run_type  UPDATE --pph_data [path/to/pph_data] \\
             --sift_run_type UPDATE --blastdb  [path/to/blastdb] \\
             --host [h] --port [p] --user [u] --pass [p] --database [db]

  General options:
    --gtf FILE           Comma-separated list of annotation GTF files; requires
                         FASTA files
    --fasta FILE         Comma-separated list of FASTA files with genomic
                         sequences; requires GTF files
    --translated FILE    Comma-separated list of FASTA files with peptide
                         sequence; skips sequence translation with GTF and FASTA
    --outdir VAL         Name of output dir (default: outdir)
    --species VAL        Latin species name (default: homo_sapiens);
                         PolyPhen-2 only works for human

  Database options (mandatory):
    --host VAL           Server host
    --port VAL           Server port
    --user VAL           Server user
    --pass VAL           Server password
    --database VAL       Name of database

  SIFT options:
    --sift_run_type VAL  SIFT run type:
                           - FULL   to run for all translations
                           - UPDATE to only run for new/changed translations
                           - NONE   to skip this step (default)
    --blastdb DIR        Path to SIFT-formatted BLAST database
                         (e.g., uniref100; required if running SIFT)
    --median_cutoff VAL  Protein alignment's median cutoff (default: 2.75)

  PolyPhen-2 options:
    --pph_run_type VAL   PolyPhen-2 run type:
                           - FULL   to run for all translations
                           - UPDATE to only run for new/changed translations
                           - NONE   to skip this step (default)
    --pph_data DIR       Path to PolyPhen-2 databases (required if running
                         PolyPhen-2); available from
                         http://genetics.bwh.harvard.edu/pph2/dokuwiki/downloads
  """
  exit 1
}

// Module imports
include { decompress }                from './nf_modules/utils.nf'
include { translate_fasta }           from './nf_modules/translations.nf'
include { store_translation_mapping } from './nf_modules/database_utils.nf'
include { run_sift_pipeline }         from './nf_modules/sift.nf'
include { run_pph2_pipeline }         from './nf_modules/polyphen2.nf'

// Check input data
if (!params.translated) {
  if (!params.fasta && !params.gtf) {
    exit 1, "ERROR: arguments --fasta/--gtf or --translated are mandatory"
  } else if (!params.fasta || !params.gtf ) {
    exit 1, "ERROR: both --fasta and --gtf need to be defined"
  }
}

if (!params.host || !params.port || !params.user || !params.pass || !params.database) {
  exit 1, "Error: --host, --port, --user, --pass and --database need to be defined"
}

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

// PolyPhen-2: check if providing databases and if species is human
if ( params.pph_run_type != "NONE" ) {
  if ( params.species != "homo_sapiens" ) {
    exit 1, "ERROR: PolyPhen-2 only works with human data"
  } else if ( !params.pph_data) {
    exit 1, "ERROR: --pph_data must be provided when running PolyPhen-2"
  }
}

// SIFT: check blastdb
if ( params.sift_run_type != "NONE" && !params.blastdb ) {
  exit 1, "ERROR: --blastdb must be supplied when running SIFT"
}

log.info """
  Predict protein function using SIFT and PolyPhen-2
  --------------------------------------------------
  species    : ${params.species}
  GTF        : ${params.gtf}
  FASTA      : ${params.fasta}
  translated : ${params.translated}
  outidr     : ${params.outdir}

  host       : ${params.host}
  port       : ${params.port}
  user       : ${params.user}
  database   : ${params.database}

  SIFT run   : ${params.sift_run_type}
  blastdb    : ${params.blastdb}

  PPH2 run   : ${params.pph_run_type}
  PPH2 data  : ${params.pph_data}
  """

def getFiles (files) {
  Channel.fromPath( files.tokenize(','), checkIfExists: true )
}

workflow {
  // Translate transcripts from GTF and FASTA if no translation FASTA is given
  if (!params.translated) {
    translate_fasta(getFiles(params.gtf), getFiles(params.fasta))
    translated = translate_fasta.out
  } else {
    translated = getFiles(params.translated)
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
