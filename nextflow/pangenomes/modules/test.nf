process test_annotation {
  // Test known annotated positions using the proper GFF3 and FASTA for the assembly
  tag "${annotation.baseName}"

  input:
    val plugin
    tuple val(accession), path(gff3), path(gff3_tbi), path(fasta), path(annotation), path(annotation_tbi)

  script:
  def output_field = plugin == 'GO' ? 'GO=' : 'PHENOTYPES='
  """
  # Build deterministic test variants from positions that are present in the
  # generated plugin annotation rather than from unrelated random exons.
  zcat ${annotation} | \
    awk -F '\t' 'BEGIN { OFS = " " }
      !/^#/ && count < 100 {
        print \$1, \$4, \$4, "C/T", "+"
        count++
      }' > input.txt
  test -s input.txt

  ${params.vep} -i input.txt \
      -o vep.out \
      --no_stats \
      --dir_plugins ${params.sw}/VEP_plugins \
      --fasta ${fasta} --gff ${gff3} \
      --plugin ${plugin},file=${annotation}

  # Confirm that VEP retrieved the expected plugin annotation.
  count=\$(awk -v key='${output_field}' '
    BEGIN { key = toupper(key) }
    !/^#/ && index(toupper(\$0), key) { count++ }
    END { print count + 0 }
  ' vep.out)
  if [ \${count} -eq 0 ]; then
    echo 'No results found with ${plugin} annotation'
    exit 1
  fi
  """
}
