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
  """
}

process run_variant_recoder {
  // Run Variant Recoder on a file with HGVS identifiers

  tag "${mappings.simpleName}"
  memory { hgvs.size() * 0.4.MB + 1.GB }

  errorStrategy 'ignore'
  maxRetries 1

  input:  tuple path(mappings), path(hgvs)
  output: tuple path(mappings), path('vr.json')

  script:
  def bin = "${params.ensembl}/ensembl-vep"
  """
  perl ${bin}/variant_recoder -i $hgvs --vcf_string > vr.json
  """
}
