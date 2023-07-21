include { fetch_licence } from '../nf_modules/fetch.nf'

workflow filter_by_licence {
  take:
    files
    licences
  main:
    fetch_licence( files )
    // Warn about discarded files
    fetch_licence.out.
      filter{ !licences.contains(it.last().text) }.
      subscribe{
        println "NOTE: Discarded ${it.first().simpleName} based on licence ${it.last().text}"
      }
    // Return files with matching licences
    filtered = fetch_licence.out.
      filter{ licences.contains(it.last().text) }.
      map{ it.first() }
  emit:
    filtered
}
