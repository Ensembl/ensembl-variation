process decompress_fasta {
  tag "${fasta.baseName}"

  input:
    tuple val(accession), path(gff3), path(fasta), path(annotation), path(annotation_tbi)
  output:
    tuple val(accession), path(gff3), path("*.fa"), path(annotation), path(annotation_tbi)

  """
  gunzip -c ${fasta} > file.fa
  """
}

process tabix_gff3 {
  tag "${gff3.baseName}"

  input:
    tuple val(accession), path(gff3), path(fasta), path(annotation), path(annotation_tbi)
  output:
    tuple val(accession), path("*.gff3.gz"), path("*.gff3.gz.tbi"), path(fasta), path(annotation), path(annotation_tbi)

  """
  gunzip -c ${gff3} | \
    awk -F '\t' 'BEGIN { OFS = "\t" }
      /^#/ { next }
      {
        if (\$9 ~ /(^|;)ID=gene:[^;]+/) {
          \$3 = "gene"
        }
        else if (\$9 ~ /(^|;)Parent=gene:[^;]+/ &&
                 \$9 ~ /(^|;)(ID=transcript:|transcript_id=)[^;]+/) {
          \$3 = "transcript"
        }

        if (\$3 != "region" && \$3 != "five_prime_UTR" && \$3 != "three_prime_UTR") {
          print
        }
      }' | \
    LC_ALL=C sort -k1,1 -k4,4n -k5,5n -t '\t' | \
    bgzip -c > file.gff3.gz
  tabix -p gff file.gff3.gz
  """
}
