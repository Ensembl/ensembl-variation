process get_hgvsp {
  // Get and sort all HGVSp identifiers from mappings file

  tag "${mappings.simpleName}"
  errorStrategy 'ignore'
  input:  path mappings
  output: tuple path(mappings), path('hgvsp.txt')

  """
  grep -Eo '".*:p..*"' $mappings |\
    sed 's/"//g' |\
    sed 's/value: //g' |\
    sort -t":" -V -k2.6 |\
    uniq > hgvsp.txt

  # error out if empty
  [[ -s hgvsp.txt ]] || exit 1
  """
}

process run_variant_recoder {
  // Run Variant Recoder on a file with HGVS identifiers

  tag "${mappings.simpleName}"
  memory { ["4 GB", "16 GB", "80 GB", "120"][task.attempt - 1] }
  maxRetries 3

  input:  tuple path(mappings), path(hgvs)
  output: tuple path(mappings), path('vr.json')

  script:
  def bin = "${ENSEMBL_ROOT_DIR}/ensembl-vep"
  """
  perl $bin/variant_recoder -i $hgvs --vcf_string > vr.json
  """
}

process parse_vr_output {
  // Prepare mappings from Variant Recoder output

  tag "${mappings.simpleName}"
  input:  tuple path(mappings), path(vr)
  output: tuple path(mappings), path('vr.txt')

  """
  parse_vr_output.py $vr -o vr.txt
  """
}
