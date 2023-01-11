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
  memory '4 GB'
  errorStrategy 'ignore'

  input:
    val peptide
    path blastdb_dir
    val blastdb_name

  output:
    path '*.alignedfasta'

  afterScript 'rm -rf *.fa *.fa.query.out'

  """
  #!/bin/csh
  cat > ${peptide.md5}.fa << EOL
${peptide.text}EOL

  setenv tmpdir "."
  setenv NCBI "/opt/blast/bin/"
  seqs_chosen_via_median_info.csh ${peptide.md5}.fa \
                                  ${blastdb_dir}/${blastdb_name} \
                                  ${params.median_cutoff}
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
  memory '4 GB'
  errorStrategy 'ignore'
  publishDir "${params.outdir}/sift"

  input:
    val peptide
    path aln

  output:
    path '*.SIFTprediction'

  afterScript 'rm -rf *.subs'

  """
  subs=${peptide.id}.subs                                                       
  create_aa_substitutions.sh sift ${peptide.id} "${peptide.seqString}" > \${subs}
  info_on_seqs ${aln} \${subs} protein.SIFTprediction
  """
}

process store_sift_scores {
  tag "${peptide.id}"
  container "ensemblorg/ensembl-vep:latest"

  input:
    val ready
    val species
    val peptide
    path weka_output

  """
  store_sift_scores.pl $species ${params.port} ${params.host} \
                       ${params.user} ${params.pass} ${params.database} \
                       ${peptide.seqString} $weka_output
  """
}

// module imports                                                               
include { delete_prediction_data; update_meta } from './database_utils.nf'        
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
  take: translated
  main:
    if ( params.sift_run_type == "UPDATE" ) {
      translated = filter_existing_translations( "sift", translated )
      wait = "ready"
    } else if ( params.sift_run_type == "FULL" ) {
      delete_prediction_data("sift")
      wait = delete_prediction_data.out
      update_sift_version()
      update_sift_db_version( file(params.blastdb) )
    }
    // Align translated sequences against BLAST database to run SIFT
    align_peptides(translated,
                   file(params.blastdb).parent,
                   file(params.blastdb).name)
    run_sift_on_all_aminoacid_substitutions(translated, align_peptides.out)
    store_sift_scores(wait, // wait for data deletion
                      params.species, translated,
                      run_sift_on_all_aminoacid_substitutions.out)
}
