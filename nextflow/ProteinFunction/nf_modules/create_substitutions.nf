#!/usr/bin/env nextflow

process getAminoacidSubstitutions {
  /*
  Generate all aminoacid substitutions for SIFT or Polyphen-2

  Returns
  -------
  Returns 1 file:
      1) Substitutions file 'subs.txt'
  */

  tag "${peptide.id}"

  input:
    val peptide
    val program 

  output:
    path 'subs.txt'

  shell:
  '''
  ALL_AAS=(A C D E F G H I K L M N P Q R S T V W Y)
  pos=0
  program=!{program}
  for ref in `grep -o . <<< !{peptide.seqString}`; do
    ((pos+=1))
    # ignore non-standard amino acids (e.g. X) when using PolyPhen-2
    if [[ $program == "polyphen2" && ! " ${ALL_AAS[*]} " =~ " ${ref} " ]]; then
      continue
    fi

    for alt in ${ALL_AAS[@]}; do
      if [[ $ref == $alt ]]; then
        continue
      elif [[ $program == "sift" ]]; then
        output="$ref$pos$alt"
      elif [[ $program == "polyphen2" ]]; then
        output="!{peptide.id}\t$pos\t$ref\t$alt"
      fi
        echo $output >> subs.txt
    done
  done
  '''
}
