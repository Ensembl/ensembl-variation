process split_by_mapping_type {
  tag "$file"

  input:
    path file

  output:
    tuple path(file), path('hgvsp.txt'), emit: hgvsp, optional: true
    tuple path(file), path('nucleotide.txt'), emit: input, optional: true

  """
  grep -m1 ":p." $file || cp $file nucleotide.txt

  if [[ ! -f nucleotide.txt ]]; then
    grep -Eo '".*:p..*"' $file | sed 's/"//g' | sed 's/value: //g' > hgvsp.txt
  fi

  # remove file if empty
  [[ -s hgvsp.txt ]] || rm -f hgvsp.txt
  """
}
