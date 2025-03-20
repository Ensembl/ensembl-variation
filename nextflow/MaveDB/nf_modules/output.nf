process concatenate_files {
  // Concatente variants associated with MaveDB scores into a single file

  input:  path(mapped_variants)
  output: path("combined.tsv")

  memory '4GB'

  """
  #!/usr/bin/env python3
  import glob, pandas, os
  from os.path import exists
  import subprocess

  output        = "combined.tsv"
  output_sorted = "combined_sorted.tsv"
  output_bgzip  = "combined.tsv.gz"
  files         = glob.glob("*map_*.tsv", recursive=True)

  print(f"Found {len(files)} files matching the pattern")

  def standardise_columns(df):
    df = df.rename(columns={'p-value':'pvalue'})
    df.columns = df.columns.str.lower()
    return df

  # concatenate header of all files
  print("Creating header...")
  header = None
  for f in files:

    content = pandas.read_csv(f, delimiter="\t", nrows=0)
    if content.empty:
      continue

    else:  
      content = standardise_columns(content)

      if header is not None:
        header = pandas.concat([header, content], axis=0, ignore_index=True)
      else:
        header = content

  print("Header columns:", header.columns.values)

  # merge data and append to file (one file at a time)
  print("\\nMerging and writing content...")
  for f in files:

    print("Processing file:", f)
    content = pandas.read_csv(f, delimiter="\t")

    if content.empty:
      continue

    else:  
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
  #!/usr/bin/env bash

  # Add hash to first line of header
  sed -i '1 s/^/#/' ${out}

  # Extract header from the combined file
  header=\$(head -n1 "${out}")

  # Remove header, LRG and chromosome patches - save in tmp.tsv
  grep -v "^#" ${out} | grep -v "^LRG" | grep -v "^CHR_" > tmp.tsv

  # Sort file by position and add header to top of file
  echo "${header}" > ${name}
  sort -k1,1 -k2,2n -k3,3n tmp.tsv | uniq >> ${name}

  # Compress and index the file
  bgzip ${name}
  tabix -s1 -b2 -e3 ${gzip}
  rm tmp.tsv
  """
}
