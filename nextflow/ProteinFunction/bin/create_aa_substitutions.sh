#!/bin/bash

program=${1?:First argument needs to be the program}
id=${2?:Second argument needs to be the ID of the peptide}
peptide=${3?:Third argument needs to be the peptide string}

pos=0
ALL_AAS=(A C D E F G H I K L M N P Q R S T V W Y)
for ref in `grep -o . <<< $peptide`; do
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
      output="$id\t$pos\t$ref\t$alt"
    else
      echo "ERROR: specified program is not supported" >&2
      exit 1
    fi
    echo -e "$output"
  done
done
