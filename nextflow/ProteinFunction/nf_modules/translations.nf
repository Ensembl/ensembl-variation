#!/usr/bin/env nextflow

/*
 * AGAT: Another GTF/GFF Analysis Toolkit
 */

process translate_fasta {
  /*
  Translate nucleotide FASTA sequences based on GTF features

  Returns
  -------
  Returns 1 file:
      1) Protein FASTA sequence 'translated.fa'
  */

  tag "${gtf} + ${fasta}"
  container "quay.io/biocontainers/agat:0.9.0--pl5321hdfd78af_0"
  label 'highmem'
  publishDir "${params.outdir}"

  input:
    path gtf
    path fasta

  output:
    path '*_translated.fa'

  """
  # decompress FASTA file if gzipped
  seq=${fasta}
  if [[ ${fasta.extension} == *gz ]]; then
    gunzip ${fasta}
    seq=${fasta.baseName}
  fi

  # decompress GTF file if gzipped
  annot=${gtf}
  if [[ ${gtf.extension} == *gz ]]; then
    gunzip ${gtf}
    annot=${gtf.baseName}
  fi

  agat_sp_extract_sequences.pl -g \${annot} -f \${seq} --protein \
                               -o ${gtf.baseName}_translated.fa
  """
}

// module imports                                                               
include { get_current_MD5_translations } from './database_utils.nf'        

workflow filter_existing_translations {
  // Filter out translation already present in database
  take:
    analysis
    translated
  main:
    get_current_MD5_translations( analysis )
    // Get MD5 hashes for translation with predictions
    current_md5s = get_current_MD5_translations.out.splitCsv().flatten().toList()
    // Concatenate MD5 hashes from both datasets; this ensures that we get MD5
    // hashes from database before running the following steps
    all_md5s = current_md5s.flatten().map{it -> [md5: it]}.concat(translated)
    // Filter out all MD5 hashes with predictions in database
    translated = all_md5s.filter { !current_md5s['value'].contains(it.md5) }
    if ( !translated.count() ) {
      exit 1, "Error: no new translations found for %analysis. All have predictions in database."
    }
  emit:
    translated
}
