process decompress_fasta {
  input:
    tuple path(gtf), path(fasta), path(annotation), path(annotation_tbi)
  output:
    tuple path(gtf), path("*.fa"), path(annotation), path(annotation_tbi)

  """
  gunzip -c ${fasta} > file.fa
  """
}

process tabix_gtf {
  input:
    tuple path(gtf), path(fasta), path(annotation), path(annotation_tbi)
  output:
    tuple path("*.gtf.gz"), path("*.gtf.gz.tbi"), path(fasta), path(annotation), path(annotation_tbi)

  """
  gunzip -c ${gtf} > file.gtf
  grep -v "#" file.gtf | sort -k1,1 -k4,4n -k5,5n -t '\t' | bgzip -c > file.gtf.gz
  tabix file.gtf.gz
  """
}
