process fetch_licence {
  // Check if files are open access (error out if not)

  errorStrategy 'ignore'

  input:  path mappings
  output: tuple path(mappings), path('LICENCE.txt')

  """
  #!/usr/bin/env python3
  import urllib.request
  import json

  urn = "urn:mavedb:${mappings.simpleName}"
  url = f"https://api.mavedb.org/api/v1/score-sets/{urn}"
  res = urllib.request.urlopen(url).read()
  res = json.loads(res)

  f = open("LICENCE.txt", "w")
  f.write(res['license']['shortName'])
  f.close()
  """
}
