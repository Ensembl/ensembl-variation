process concatenate_files {
  // Concatente variants associated with MaveDB scores into a single file

  input:  path(mapped_variants)
  output: path("combined.tsv")
  """
  #!/usr/bin/env python3
  import glob, pandas
  from os.path import exists
  import subprocess

  output        = "combined.tsv"
  output_sorted = "combined_sorted.tsv"
  output_bgzip  = "combined.tsv.gz"

  # concatenate header of all files
  print("Creating header...")
  header = None
  for f in glob.glob("*.tsv"):
    content = pandas.read_csv(f, delimiter="\t", nrows=0)
    if header is not None:
      header = pandas.concat([header, content], axis=0, ignore_index=True)
    else:
      header = content
  print(header)

  # merge data and append to file (one file at a time)
  print("Merging and writing content...")
  for f in glob.glob("*.tsv"):
    print(f)
    content = pandas.read_csv(f, delimiter="\t")
    out = pandas.concat([header, content], axis=0, ignore_index=True)
    out.to_csv(output, sep="\t", mode="a", index=False, header=not exists(output))
  """
}

process tabix {
  input:
    path out
  output:
    path "${name}"
    path "${name}.tbi"

  //publishDir ${params.output}, mode: 'move', overwrite: true

  script:
  def name="MaveDB_variants.tsv"
  def gzip=name + ".gz"
  """
  sort -k1 -nk2,3 ${out} > ${name}
  bgzip ${name}
  tabix -s1 -b2 -e3 ${gzip}

  rm ${name}
  """
}
