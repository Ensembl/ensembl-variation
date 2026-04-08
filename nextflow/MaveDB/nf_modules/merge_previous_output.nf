process merge_previous_output {
  tag "merge_previous_output"

  input:
    tuple path(current_combined), path(prev_output), path(skipped_urns)

  output:
    path("merged_combined.tsv")

  script:
  """
  #!/usr/bin/env python3
  import csv
  import gzip

  current_combined = "${current_combined}"
  prev_output = "${prev_output}"
  skipped_urns = "${skipped_urns}"
  merged_output = "merged_combined.tsv"

  def normalise_header(columns):
      normalised = []
      for idx, col in enumerate(columns):
          value = col.strip()
          if idx == 0 and value.startswith("#"):
              value = value[1:]
          value = value.lower()
          if value == "p-value":
              value = "pvalue"
          normalised.append(value)
      return normalised

  def read_tsv(path, opener):
      with opener(path, "rt", newline="") as handle:
          reader = csv.reader(handle, delimiter="\\t")
          rows = list(reader)
      if not rows:
          raise SystemExit(f"ERROR: empty TSV input: {path}")
      header = normalise_header(rows[0])
      return header, rows[1:]

  def project_rows(header, rows, union_header):
      index = {name: pos for pos, name in enumerate(header)}
      projected = []
      for row in rows:
          projected.append(tuple(row[index[name]] if name in index and index[name] < len(row) else "" for name in union_header))
      return projected

  if not skipped_urns:
      raise SystemExit("ERROR: skipped URN file path not provided")

  with open(skipped_urns) as handle:
      skipped = {line.strip() for line in handle if line.strip()}

  current_header, current_rows = read_tsv(current_combined, open)
  if "urn" not in current_header:
      raise SystemExit(f"ERROR: urn column not found in current file: {current_combined}")

  previous_header, previous_rows = read_tsv(prev_output, gzip.open)
  if "urn" not in previous_header:
      raise SystemExit(f"ERROR: urn column not found in previous output: {prev_output}")

  urn_idx_prev = previous_header.index("urn")
  filtered_previous_rows = [row for row in previous_rows if urn_idx_prev < len(row) and row[urn_idx_prev] in skipped]

  union_header = list(current_header)
  for column in previous_header:
      if column not in union_header:
          union_header.append(column)

  merged_rows = []
  seen = set()
  for row in project_rows(current_header, current_rows, union_header) + project_rows(previous_header, filtered_previous_rows, union_header):
      if row not in seen:
          seen.add(row)
          merged_rows.append(row)

  with open(merged_output, "w", newline="") as handle:
      writer = csv.writer(handle, delimiter="\\t", lineterminator="\\n")
      writer.writerow(union_header)
      writer.writerows(merged_rows)
  """
}
