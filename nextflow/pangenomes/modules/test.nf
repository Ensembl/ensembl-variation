process test_annotation {
  // Test annotation using the proper GTF and FASTA for the assembly
  tag "${annotation.baseName}"

  input:
    val plugin
    tuple path(gtf), path(gtf_tbi), path(fasta), path(annotation), path(annotation_tbi)

  """
  # create random examples based on GTF exons
  zgrep exon ${gtf} | \\
    awk '\$3 == "exon"' | \\
    head -n 10000 | \\
    shuf -n 100 | \\
    awk '{print \$1, \$4, \$4, "C/T", "+"}' > input.txt
  
  ${params.vep} -i input.txt \\
      -o vep.out \\
      --fasta ${fasta} --gtf ${gtf} \\
      --plugin ${plugin},file=${annotation}

  # count number of lines with plugin annotation
  count=\$(grep -v "^#" vep.out | grep -c -i ${plugin})
  if [ \${count} -eq 0 ]; then
    echo 'No results found with ${plugin} annotation'
    exit 1
  fi
  """
}
