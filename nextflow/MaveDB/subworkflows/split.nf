process inspect_mappings {
  tag "${urn}"
  errorStrategy 'terminate'

  input:
    tuple val(urn), path(mappings), path(scores), path(metadata)
  output:
    tuple val(urn), path(mappings), path(scores), path(metadata), path('mapping_type.txt'), path('hgvsp.txt')

  """
  mapping_hgvs.py inspect --hgvsp-output hgvsp.txt "${mappings}" > mapping_type.txt
  """
}

workflow split_by_mapping_type {
  take:
    files

  main:
    inspected = inspect_mappings(files.map {
      tuple(it.urn, it.mappings, it.scores, it.metadata)
    })
    type = inspected.branch {
      hgvs_pro: it[4].text.trim() == "hgvs.p"
      hgvs_nt:  it[4].text.trim() == "hgvs.g"
      unmapped: it[4].text.trim() == "unmapped"
    }

  emit:
    hgvs_pro = type.hgvs_pro.map { urn, mappings, scores, metadata, mapping_type, hgvsp ->
      tuple(urn, mappings, scores, metadata, hgvsp)
    }
    hgvs_nt = type.hgvs_nt.map { urn, mappings, scores, metadata, mapping_type, hgvsp ->
      tuple(urn, mappings, scores, metadata, mapping_type.text.trim())
    }
    unmapped = type.unmapped.map { urn, mappings, scores, metadata, mapping_type, hgvsp ->
      tuple(urn, mappings, scores, metadata, "unmapped")
    }
}
