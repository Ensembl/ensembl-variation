process download_MaveDB_metadata {
  // Download MaveDB metadata
  errorStrategy 'ignore'

  input:  val urn
  output: tuple val(urn), path('metadata.json'), path('LICENCE.txt')

  """
  #!/usr/bin/env python3
  import urllib.request
  import json

  url = "https://api.mavedb.org/api/v1/score-sets/${urn}"
  res = urllib.request.urlopen(url).read()
  res = json.loads(res)
  
  with open('metadata.json', 'w') as f:
    json.dump(res, f)

  with open("LICENCE.txt", "w") as f:
    f.write(res['license']['shortName'])
  """
}

process download_MaveDB_data {
  // Download MaveDB data
  input: tuple val(urn), path(metadata), path(licence)
  output: tuple val(urn), path('mappings.json'), path('scores.csv'), path(metadata)

  script:
    url = "https://api.mavedb.org/api/v1/score-sets/${urn}"
  """
  wget ${url}/mapped-variants -O mappings.json
  wget ${url}/scores -O scores.csv
  """
}