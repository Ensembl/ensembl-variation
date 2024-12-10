#!/usr/bin/env nextflow

/* 
 * Predict protein function using SIFT
 */

process get_sift_version {
  container "ensemblorg/sift:6.2.1"
  output: stdout
  """
  grep -Po "VERSION \\d+\\.\\d+\\.\\d+" /opt/sift/VERSION_UPDATE | tail -n 1 | tr -d '\n'
  """
}

process get_sift_db_version {
  /*
  Get version from UniRef release note file (if available)
  */
  container "ensemblorg/sift:6.2.1"
  input:
    path db_dir
    path db
  output:
    stdout

  """
  release_note=${db_dir}/${db.baseName}.release_note
  if [ -f \${release_note} ]; then
    version=`grep -oE "[0-9]+_[0-9]+" \${release_note}`
  else
    version=`date -r ${db} -u +%Y_%m`
  fi
  echo -n "${db.baseName} (\${version})"
  """
}

process align_peptides {
  /*
  Run multiple alignment for peptide sequence
  */

  tag "${peptide.md5}"
  container "ensemblorg/sift:6.2.1"

  memory { peptide.size() * 40.MB + 4.GB }
  errorStrategy 'ignore'
  maxRetries 1

  input:
    val peptide
    path blastdb_dir
    val blastdb_name

  output:
    tuple val(peptide), path('*.alignedfasta'), emit: aln optional true
    path 'expected_error.txt', emit: errors optional true

  afterScript 'rm -rf *.fa *.fa.query.out'

  """
  #!/bin/csh
  cat > ${peptide.md5}.fa << EOL
${peptide.text}EOL

  set error=0
  setenv tmpdir "."
  #setenv NCBI "/opt/blast/bin/"
  setenv NCBI "/hps/software/users/ensembl/ensw/C8-MAR21-sandybridge/linuxbrew/Cellar/blast/2.2.30/bin"
  seqs_chosen_via_median_info.csh ${peptide.md5}.fa \\
                                  ${blastdb_dir}/${blastdb_name} \\
                                  ${params.median_cutoff} || setenv error \$status
  # Capture expected errors to avoid failing job
  set nonomatch=1
  if ( -f *.error ) capture_expected_errors.sh ${peptide.md5} sift_align *.error \$error expected_error.txt
  """
}

process run_sift_on_all_aminoacid_substitutions {
  /*
  Run SIFT

  Returns
  -------
  Returns 1 file:
      1) Output 'protein.SIFTprediction'
  */

  tag "${peptide.md5}"
  container "ensemblorg/sift:6.2.1"
  publishDir "${params.outdir}/sift"

  memory { peptide.size() * 40.MB + 4.GB }
  errorStrategy 'ignore'
  maxRetries 1

  input:
    tuple val(peptide), path(aln)

  output:
    tuple val(peptide), path('*.SIFTprediction'), emit: scores optional true
    path 'expected_error.txt', emit: errors optional true

  afterScript 'rm -rf *.subs'

  """
  subs=${peptide.id}.subs                                                       
  create_aa_substitutions.sh sift ${peptide.id} "${peptide.seqString}" > \${subs}

  error=0
  info_on_seqs ${aln} \${subs} protein.SIFTprediction || error=\$?

  # Capture expected errors to avoid failing job
  capture_expected_errors.sh ${peptide.md5} sift .command.err \$error expected_error.txt
  """
}

process store_sift_scores {
  tag "${peptide.id}"
  container "ensemblorg/ensembl-vep:latest"
  time '10m'

  cache false

  input:
    val ready
    val species
    tuple val(peptide), path(sift_scores)

  """
  store_sift_scores.pl ${species} ${params.offline} ${params.sqlite_db} \
                       ${params.port} ${params.host} ${params.user} ${params.pass} ${params.database} \
                       ${peptide.seqString} ${sift_scores}
  """
}

// module imports                                                               
include { delete_prediction_data; update_meta } from './database.nf'        
include { filter_existing_translations        } from './translations.nf'

workflow update_sift_version {
  get_sift_version()
  update_meta("sift_version", get_sift_version.out)
}

workflow update_sift_db_version {
  take: db
  main:
    get_sift_db_version( db.parent, db )
    update_meta("sift_db_version", get_sift_db_version.out)
}

workflow run_sift_pipeline {
  take: 
    translated
    sqlite_db_prep
  main:
    if ( params.sift_run_type == "UPDATE" && !params.offline ) {
      translated = filter_existing_translations( "sift", translated )
      wait = "ready"
    } else if ( params.sift_run_type == "FULL" && !params.offline ) {
      delete_prediction_data("sift")
      wait = delete_prediction_data.out
      update_sift_version()
      update_sift_db_version( file(params.blastdb) )
    }
    // Align translated sequences against BLAST database to run SIFT
    blast = align_peptides(translated,
                           file(params.blastdb).parent,
                           file(params.blastdb).name)
    sift = run_sift_on_all_aminoacid_substitutions(blast.aln)
    store_sift_scores(wait, // wait for data deletion
                      params.species, sift.scores)
    all_errors = sift.errors.concat(blast.errors)
  emit:
    errors = all_errors
}
