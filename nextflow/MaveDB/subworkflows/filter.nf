include { download_MaveDB_metadata } from '../nf_modules/fetch.nf'

workflow filter_by_licence {
  take:
    urn
    licences
  main:
    download_MaveDB_metadata( urn )
    // Warn about discarded files
    data = download_MaveDB_metadata.out
             .map { [ urn: it[0], metadata: it[1], licence: it[2] ] }

    data.
      filter { !licences.contains(it.licence.text) }.
      subscribe {
        println "NOTE: Discarded ${it.urn} based on licence ${it.licence.text}"
      }

    // Return files with matching licences
    filtered = data.filter { licences.contains(it.licence.text) }
  emit:
    filtered
}
