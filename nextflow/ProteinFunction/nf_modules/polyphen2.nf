#!/usr/bin/env nextflow

/* 
 * Predict protein function using PolyPhen
 */

params.pph_data = "/hps/nobackup/flicek/ensembl/variation/nuno/sift-polyphen2-nextflow-4667/input/polyphen2"

process get_pph2_version {
  container "nunoagostinho/polyphen-2:2.2.3"
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
  container "nunoagostinho/polyphen-2:2.2.3"
  containerOptions "--bind ${params.pph_data}:/opt/pph2/data"
  memory '4 GB'
  errorStrategy 'ignore'

  input:
    val peptide

  output:
    path '*_scores.txt'

  shell:
  '''
  subs=!{peptide.id}.subs
  create_aa_substitutions.sh polyphen2 !{peptide.id} \
                             "!{peptide.seqString}" > ${subs}

  fasta=!{peptide.id}.fa
  cat > ${fasta} <<EOL
!{peptide.text}EOL

  mkdir -p tmp/lock
  out=!{peptide.id}_scores.txt
  /opt/pph2/bin/run_pph.pl -A -d tmp -s ${fasta} ${subs} 1> $out

  # Remove output if only contains header
  if [ "$( wc -l <$out )" -eq 1 ]; then rm %out; fi
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

  //tag "${pph2_out.baseName} $model"
  container "nunoagostinho/polyphen-2:2.2.3"
  errorStrategy 'ignore'

  input:
    each model
    path pph2_out

  output:
    path '*.txt', emit: results
    path '*.out', emit: processed
    path '*.err', emit: error
    val "${model}", emit: model

  shell:
  '''
  res=!{pph2_out.baseName}_!{model}.txt
  run_weka.pl -l /opt/pph2/models/!{model} !{pph2_out} \
              1> $res 2> !{pph2_out.baseName}_!{model}.err

  # clean results and append variant coordinates and transcript id
  var=$( echo !{pph2_out.baseName}-PolyPhen2 | sed 's/-/,/g' )
  grep -v "^#" ${res} | awk -v var="$var" -v OFS=',' -F'\t' '{$1=$1;print var,$12,$16}' | sed 's/,/\t/g' > !{pph2_out.baseName}.out
  '''
}

process store_pph2_scores {
  tag "${peptide.id}"
  container "ensemblorg/ensembl-vep:latest"
  input:
    val ready
    val species
    val peptide
    path weka_output
    val model

  """
  store_polyphen_scores.pl $species ${params.port} ${params.host} \
                           ${params.user} ${params.pass} ${params.database} \
                           ${peptide.seqString} $weka_output $model
  """
}

// module imports                                                               
include { delete_prediction_data; update_meta } from './database_utils.nf'        
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
  run_pph2_on_all_aminoacid_substitutions(translated)

  weka_model = Channel.of("HumDiv.UniRef100.NBd.f11.model",
                          "HumVar.UniRef100.NBd.f11.model")
  run_weka(weka_model, run_pph2_on_all_aminoacid_substitutions.out)
  store_pph2_scores(wait, // wait for data deletion
                    params.species, translated,
                    run_weka.out.results, run_weka.out.model)
}
