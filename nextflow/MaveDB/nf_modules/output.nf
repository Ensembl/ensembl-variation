process concatenate_files {
  // Concatente variants associated with MaveDB scores into a single file

  input:  path(mapped_variants)
  output: path("combined.tsv")

  memory '4GB'

  """
  #!/usr/bin/env python3
  import glob, pandas
  from os.path import exists
  import subprocess

  output        = "combined.tsv"
  output_sorted = "combined_sorted.tsv"
  output_bgzip  = "combined.tsv.gz"
  files         = "*.tsv"

  def standardise_columns(df):
    df = df.rename(columns={'p-value':'pvalue'})
    df.columns = df.columns.str.lower()
    return df

  # concatenate header of all files
  print("Creating header...")
  header = None
  for f in glob.glob(files):
    content = pandas.read_csv(f, delimiter="\t", nrows=0)
    content = standardise_columns(content)

    if header is not None:
      header = pandas.concat([header, content], axis=0, ignore_index=True)
    else:
      header = content
  print(header.columns.values)

  # merge data and append to file (one file at a time)
  print("\\nMerging and writing content...")
  for f in glob.glob(files):
    print(f)
    content = pandas.read_csv(f, delimiter="\t")
    content = standardise_columns(content)
    out = pandas.concat([header, content], axis=0, ignore_index=True)
    out.to_csv(output, sep="\t", mode="a", index=False, header=not exists(output))
  """
}

process tabix {
  publishDir file(params.output).parent, mode: 'move', overwrite: true

  input:  path out
  output: path "*"

  script:
  def name = file(params.output).baseName
  def gzip = file(params.output).name
  """
  # add hash to first line of header
  sed -i '1 s/^/#/' ${out}

  # remove LRG and chromosome patches
  grep -v "^#" ${out} | grep -v "^LRG" | grep -v "^CHR_" > tmp.tsv

  # sort file by position
  (head -n1 ${out}; sort -k1,1 -k2,2n -k3,3n tmp.tsv | uniq) > ${name}

  bgzip ${name}
  tabix -s1 -b2 -e3 ${gzip}
  rm tmp.tsv
  """
}
