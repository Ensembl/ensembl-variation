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
