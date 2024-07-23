process run_variant_recoder {
  // Run Variant Recoder on a file with HGVS identifiers
  label 'bigmem'

  input:  tuple val(urn), path(mappings), path(scores), path(metadata), path(hgvs)
  output: tuple val(urn), path(mappings), path(scores), path(metadata), path('vr.json')

  tag "${urn}"
  memory { file(hgvs.target).countLines() * 100.MB + 4.GB }

  errorStrategy 'ignore'
  maxRetries 1

  script:
  def bin = "${params.ensembl}/ensembl-vep"
  def reg = params.registry ? "--registry ${params.registry}" : ""
  """
  perl ${bin}/variant_recoder -i $hgvs --vcf_string $reg > vr.json
  """
}
