process get_licence {
  // Check if files are open access (error out if not)

  errorStrategy 'ignore'

  input:  path mappings
  output: tuple path(mappings), path('LICENCE.txt')

  """
  #!/usr/bin/env python3
  import urllib.request
  import json

  urn = "urn:mavedb:${mappings.simpleName}"
  url = f"https://api.mavedb.org/api/v1/scoresets/{urn}"
  res = urllib.request.urlopen(url).read()
  res = json.loads(res)

  f = open("LICENCE.txt", "w")
  f.write(res['license']['shortName'])
  f.close()
  """
}

workflow get_files_by_licence {
  take:
    files
    licences
  main:
    get_licence( files )
    // Warn about discarded files
    get_licence.out.
      filter{ !licences.contains(it.last().text) }.
      subscribe{
        println "NOTE: Discarded ${it.first().simpleName} based on licence ${it.last().text}"
      }
    // Return files with matching licences
    filtered = get_licence.out.
      filter{ licences.contains(it.last().text) }.
      map{ it.first() }
  emit:
    filtered
}

def split_by_mapping_type (files) {
  // split mapping files based on HGVS type (HGVSp or HGVSg files)
  type = files.map {
    it.withReader {
      while( line = it.readLine() ) {
        if (line.contains("hgvs.")) {
          // get first line describing HGVS type
          if (line.contains("hgvs.p")) {
            hgvs = "hgvs.p"
          } else if (line.contains("hgvs.g")) {
            hgvs = "hgvs.g"
          } else {
            throw new Exception("Error: HGVS type in '${line.trim()}' not expected")
          }
          break
        }
      }
    }
    [file: it, hgvs: hgvs]
  }.branch{
    hgvs_pro: it.hgvs == "hgvs.p"
    hgvs_nt:  it.hgvs == "hgvs.g"
  }

  // clean up: only return the files in each branch
  files = [hgvs_pro: null, hgvs_nt: null]
  files.hgvs_pro = type.hgvs_pro.map { it.file }
  files.hgvs_nt  = type.hgvs_nt.map { it.file }
  return files
}
