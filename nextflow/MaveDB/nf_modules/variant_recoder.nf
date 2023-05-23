process get_hgvsp {
  // Get and sort all HGVSp identifiers from mappings file

  tag "$mappings"
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

  tag "$file"
  input:  tuple path(file), path(hgvs)
  output: tuple path(file), path('vr.json')

  """
  variant_recoder -i $hgvs --vcf_string > vr.json
  """
}

process prepare_vr_mappings {
  // Prepare mappings from Variant Recoder output

  tag "$file"
  input:  tuple path(file), path(vr)
  output: tuple path(file), path('vr.txt')

  """
  #!/usr/bin/env python3
  import csv
  import json

  f = open('$vr')
  data = json.load(f)

  line = []
  for result in data:
    for allele in result:
      info = result[allele]

      if type(info) is list and "Unable to parse" in info[0]:
        continue

      hgvsp  = info["input"]
      for string in info["vcf_string"]:
        chr, pos, ref, alt = string.split('-')
        line.append( [hgvsp, chr, pos, ref, alt] )
  f.close()

  w = open('vr.txt', 'w')
  writer = csv.writer(w, delimiter='\t')
  writer.writerows(line)
  w.close()
  """
}
