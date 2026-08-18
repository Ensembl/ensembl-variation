def split_by_mapping_type (files) {
  // split mapping files based on HGVS type (HGVSp or HGVSg files)
  type = files.map {
    def hgvs = "unmapped"
    def stdout = new StringBuffer()
    def stderr = new StringBuffer()
    def proc = ["${baseDir}/bin/mapping_hgvs.py", "type", it.mappings.toString()].execute()
    proc.waitForProcessOutput(stdout, stderr)

    if (proc.exitValue() == 0) {
      hgvs = stdout.toString().trim()
    } else {
      println "No current mapped HGVS expression found for ${it.urn}: ${stderr.toString().trim()}"
    }

    it + [hgvs: hgvs]
  }.branch{
    hgvs_pro: it.hgvs == "hgvs.p"
    hgvs_nt:  it.hgvs == "hgvs.g"
    unmapped: it.hgvs == "unmapped"
  }
  return type
}
