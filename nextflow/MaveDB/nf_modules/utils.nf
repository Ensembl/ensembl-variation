process get_hgvsp {
  // Get and sort all HGVSp identifiers from mappings file

  tag "${urn}"
  errorStrategy 'ignore'
  input:  tuple val(urn), path(mappings), path(scores), path(metadata), val(hgvs)
  output: tuple val(urn), path(mappings), path(scores), path(metadata), path('hgvsp.txt')

  """
  sed 's/"/\\n/g' $mappings |\\
    grep -Eo '.*:p\\..*' |\\
    sort -t":" -V -k2.6 |\\
    uniq > hgvsp.txt
  """
}
