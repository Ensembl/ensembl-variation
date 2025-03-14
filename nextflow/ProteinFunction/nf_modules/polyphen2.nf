#!/usr/bin/env nextflow

/* 
 * Predict protein function using PolyPhen
 */

params.pph_data = "/hps/nobackup/flicek/ensembl/variation/nuno/sift-polyphen2-nextflow-4667/input/polyphen2"

process get_pph2_version {
  container "ensemblorg/polyphen-2:2.2.3"
  output: stdout

  """
  cat /opt/pph2/VERSION | tr -d '\n'                                           
  """
}

process run_pph2_on_all_aminoacid_substitutions {
  /*
  Run PolyPhen-2 on a protein sequence with a substitions file

  Returns
  -------
  Returns 2 files:
      1) Output '*.txt'
      2) Errors '*.err'
  */

  tag "${peptide.md5}"
  container "ensemblorg/polyphen-2:2.2.3"
  containerOptions "--bind ${params.pph_data}:/opt/pph2/data"

  memory { peptide.size() * 100.MB + 4.GB }
  errorStrategy 'ignore'
  maxRetries 1

  input:
    val peptide

  output:
    tuple val(peptide), path ('*_scores.txt'), emit: scores optional true
    path 'expected_error.txt', emit: errors optional true

  afterScript 'rm -rf *.fa *.subs tmp/'

  shell:
  '''
  subs=!{peptide.id}.subs
  create_aa_substitutions.sh polyphen2 !{peptide.id} \
                             "!{peptide.seqString}" > ${subs}

  fasta=!{peptide.id}.fa
  cat > ${fasta} <<EOL
!{peptide.text}EOL

  error=0
  mkdir -p tmp/lock
  out=!{peptide.id}_scores.txt
  /opt/pph2/bin/run_pph.pl -A -d tmp -s ${fasta} ${subs} > $out || error=$?

  # Capture expected errors to avoid failing job
  capture_expected_errors.sh !{peptide.md5} pph2 .command.err $error expected_error.txt

  # Remove output if only contains header
  if [ "$( wc -l <$out )" -eq 1 ]; then rm $out; fi
  '''
}

process run_weka {
  /*
  Run Weka

  Returns
  -------
  Returns 2 files:
      1) Output '*.txt'
      2) Error '*.err'
  */

  tag "${peptide.md5} ${model}"
  container "ensemblorg/polyphen-2:2.2.3"
  label 'retry_before_ignoring'

  input:
    each model
    tuple val(peptide), path(pph2_out)

  output:
    tuple val(peptide), path('*.txt'), val("${model}")

  """
  run_weka.pl -l /opt/pph2/models/${model} ${pph2_out} \
              > ${pph2_out.baseName}_${model}.txt
  """
}

process store_pph2_scores {
  tag "${peptide.md5} ${model}"
  container "ensemblorg/ensembl-vep:latest"
  time '10m'

  cache false

  input:
    val ready
    val species
    tuple val(peptide), path(weka_output), val(model)

  """
  store_polyphen_scores.pl $species ${params.port} ${params.host} \
                           ${params.user} ${params.pass} ${params.database} \
                           ${peptide.seqString} ${weka_output} ${model}
  """
}

// module imports                                                               
include { delete_prediction_data; update_meta } from './database.nf'        
include { filter_existing_translations        } from './translations.nf'

workflow run_pph2_pipeline {
  take: translated
  main:
    if ( params.pph_run_type == "UPDATE" ) {
      translated = filter_existing_translations( "polyphen_%", translated )
      wait = "ready"
    } else if ( params.pph_run_type == "FULL" ) {
      delete_prediction_data("polyphen_%")
      wait = delete_prediction_data.out
      get_pph2_version()
      update_meta("polyphen_version", get_pph2_version.out)
    }
    // Run PolyPhen-2 and Weka
    pph2 = run_pph2_on_all_aminoacid_substitutions(translated)

    weka_model = Channel.of("HumDiv.UniRef100.NBd.f11.model",
                            "HumVar.UniRef100.NBd.f11.model")
    weka = run_weka(weka_model, pph2.scores)
    store_pph2_scores(wait, // wait for data deletion
                      params.species, weka)
  emit:
    errors = pph2.errors
}
