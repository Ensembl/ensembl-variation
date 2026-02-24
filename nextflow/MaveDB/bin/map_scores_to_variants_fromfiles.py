#!/usr/bin/env python3
import urllib.request                     # For making HTTP requests (e.g. to the Ensembl REST API)
from urllib.parse import urlencode        # To URL‐encode query parameters
import csv                                # For reading/writing CSV/TSV files
import os
import json                               # For parsing JSON files (mappings, metadata, etc.)
import operator                           # (Imported but not used explicitly)
import warnings                           # For issuing warnings with custom formatting
from collections import OrderedDict       # To maintain key order in output dictionaries
from math import isclose                  # (Imported but not used explicitly)
import re                                 # (Imported but not used explicitly)
import argparse                           # For parsing command-line arguments
import sys                                # For writing to stderr
import datetime                           # For timestamping structured logs

# Structured logging (keeps stdout clean and uniform across all tools)
URN  = os.environ.get("MAVEDB_URN", "na")
STEP = os.environ.get("STEP", "map_scores")
def log(reason, subid="na", **kv):
    """
    Uniform structured stderr logger to match the pipeline’s format.
    """
    ts = datetime.datetime.now().astimezone().isoformat()
    extra = " ".join(f"{k}={json.dumps(v, ensure_ascii=False)}" for k, v in kv.items())
    sys.stderr.write(f"[{ts}][MaveDB][URN={URN}][STEP={STEP}][REASON={reason}][SUBID={subid}] {extra}\n")
    sys.stderr.flush()

# Customize warning messages: route Python warnings through the structured logger
def customshowwarning(message, category, filename, lineno, file=None, line=None):
    log("warning", msg=str(message))
warnings.showwarning = customshowwarning

def main():
  # Set up argument parsing for input files and parameters
  parser = argparse.ArgumentParser(
    description='Output file with MaveDB scores mapped to variants')
  parser.add_argument('--vr', type=str,
                      help="path to file containg Variant Recoder output with 'vcf_string' enabled (optional)")
  parser.add_argument('--urn', type=str, help="MaveDB URN")
  parser.add_argument('--scores', type=str,
                      help="path to file with MaveDB URN scores")
  parser.add_argument('--mappings', type=str,
                      help="path to file with MaveDB URN mappings")
  parser.add_argument('--metadata', type=str,
                      help="path to file with MaveDB URN metadata")
  parser.add_argument('-o', '--output', type=str,
                      help="path to output file")
  parser.add_argument('--round', type=int,
                      help="Number of decimal places for rounding values (default: not used)")
  args = parser.parse_args()

  # Load the scores (CSV), mappings (JSON), and metadata (JSON) files
  log("load_start", scores=args.scores, mappings=args.mappings, metadata=args.metadata, vr=args.vr)

  scores = load_scores(args.scores)

  # Check if the scores file is empty, if so, log an error and exit - don't process this URN
  if not scores:
    log("scores_empty", target=args.scores, urn=args.urn)
    sys.exit(1)

  # Throw an info/warning if scores file has very few entries, as this may explain lack of mappings
  if len(scores) < 50:  # (original had 10; using 50 to surface small datasets earlier)
    log("scores_small", n=len(scores), urn=args.urn)

  with open(args.mappings) as f:
    mappings = json.load(f)

  with open(args.metadata) as f:
    metadata = json.load(f)

  # Extract MaveDB IDs from each mapping entry (robust access)
  mapped_scores = mappings.get("mapped_scores", [])
  mavedb_ids = [i.get("mavedb_id") for i in mapped_scores if "mavedb_id" in i]

  # Pre-build a dictionary mapping accession IDs to mapping records
  # This makes lookups robust and order-independent
  mapping_dict = {m["mavedb_id"]: m for m in mapped_scores if "mavedb_id" in m}

  log("mappings_loaded", n=len(mapped_scores), unique_ids=len(mapping_dict))

  # If a Variant Recoder output file is provided, load it; otherwise, set matches to None
  if args.vr is not None:
    hgvsp2vars = load_vr_output(args.vr)
  else:
    hgvsp2vars = None

  # Create the mapping between variant coordinates and MaveDB scores
  log("map_prepare")
  mapped_data = map_scores_to_variants(scores, mappings, metadata, mavedb_ids, hgvsp2vars, args.round, mapping_dict)
  write_variant_mapping(args.output, mapped_data)

  log("done", out=args.output)
  return True

def load_vr_output (f):
  """
  Load Variant Recoder output.

  Iterates over each allele in the JSON data, skipping any entries that couldn't be parsed.
  For each valid allele, splits each VCF string (format: chr-start-ref-alt) and builds an
  ordered dictionary with variant details. Returns a dictionary mapping HGVS strings to lists of variants.
  """
  data = json.load(open(f))

  # Check if file has warnings - may indicate that variant recoder couldn't recode some or all variants
  if any("warnings" in item for item in data):
    log("vr_has_warnings", file=f)
    # If *all* entries are warnings => nothing parseable
    if all("warnings" in item for item in data):
      log("vr_only_unable_to_parse", file=f)
      sys.exit(1)

  matches = {}
  for result in data:
    for allele in result:
      info = result[allele]
      # Skip and log entries that variant_recoder could not parse
      if isinstance(info, list) and (info and ("Unable to parse" in info[0] or "skipped" in info[0])):
        log("vr_unable_to_parse", subid=allele, msg=info[0] if info else "")
        continue
      hgvs = info["input"]
      for string in info["vcf_string"]:
        chr, start, ref, alt = string.split('-')
        end = int(start) + len(alt) - 1
        d = OrderedDict([("HGVSp", hgvs),
                         ("chr",   chr),
                         ("start", start),
                         ("end",   end),
                         ("ref",   ref),
                         ("alt",   alt)])
        if hgvs not in matches:
          matches[hgvs] = []
        matches[hgvs].append(d)
  log("vr_loaded", n=sum(len(v) for v in matches.values()))
  return matches

def load_HGVSp_to_variant_matches (f):
  """Load HGVSp to variant matches from a TSV file."""
  matches = {}
  with open(f) as csvfile:
    reader = csv.DictReader(csvfile, delimiter="\t")
    for row in reader:
      hgvs = row['HGVSp']
      if hgvs not in matches:
        matches[hgvs] = []
      matches[hgvs].append(row)
  return(matches)

def load_scores (f):
  """Load MaveDB scores from a CSV file into a list of dictionaries."""
  scores = []
  with open(f) as csvfile:
    reader = csv.DictReader(csvfile)
    # Strip whitespace from each header name -- I think only necessary due to the csv viewer adding spacing and then this was cached in a nf run. Consider removing.
    reader.fieldnames = [field.strip() for field in reader.fieldnames]
    for row in reader:
      # Strip whitespace from each value if it is a string -- same note as above
      clean_row = { key: value.strip() if isinstance(value, str) else value for key, value in row.items() }
      scores.append(clean_row)
  log("scores_loaded", n=len(scores))
  return scores

# Global variable for caching chromosome name
chrom = None
def get_chromosome (hgvs):
  """
  Lookup chromosome name in the Ensembl REST API (caches result globally).

  Splits the HGVS string to extract the chromosome, retrieves synonyms from Ensembl,
  and selects the UCSC name (removing 'chr' if present).
  """
  global chrom
  if (chrom is None):
    chrom = hgvs.split(":")[0]
    url   = f"https://rest.ensembl.org/info/assembly/homo_sapiens/{chrom}?"
    data  = urlencode({"synonyms": 1, "content-type": "application/json"})
    res   = urllib.request.urlopen(url + data).read()
    res   = json.loads(res)
    chrom = [each["name"] for each in res["synonyms"] if each['dbname'] == "UCSC"][0]
    chrom = chrom.replace("chr", "")
  return chrom

def join_information(hgvs, mapped_info, row, extra):
  """
  Join variant and MaveDB score information for a given HGVS.

  It extracts:
    - Start and end coordinates from 'location'.
    - The reference allele from the first element of 'extensions'.
    - The alternate allele from 'state'.
    - The HGVS expression from the first element of 'expressions'.

  Then, it builds and returns an ordered dictionary containing these values merged with extra metadata and the score row.
  """
  # Extract coordinate information.
  start = mapped_info["location"]["start"]
  end   = mapped_info["location"]["end"]
  # Extract the reference allele.
  ref   = mapped_info["extensions"][0]["value"]
  # Extract the alternate allele.
  alt   = mapped_info["state"]["sequence"]
  # Extract the HGVS expression.
  hgvs  = mapped_info["expressions"][0]["value"]

  mapped = OrderedDict([
    ("chr", get_chromosome(hgvs)),  # Lookup the chromosome using the full HGVS string.
    ("start", start + 1),           # Convert 0-based to 1-based indexing.
    ("end", end),
    ("ref", ref),
    ("alt", alt),
    ("hgvs", hgvs)
  ])

  # Merge extra metadata and the current score row.
  mapped.update(extra)
  mapped.update(row)
  return [mapped]

def match_information (hgvs, matches, row, extra):
  """
  Match a given HGVS to variant details using pre-loaded HGVSp-variant matches.

  If the provided HGVS is not found in the matches, a log entry is emitted.
  For each matching entry, merge the match with extra metadata and the score row.
  """
  out = []
  if hgvs not in matches:
    log("hgvsp_not_in_vr_matches", subid=hgvs)
    return out

  for match in matches[hgvs]:
    mapped = match
    mapped['hgvs'] = match['HGVSp']
    mapped.update(extra)
    mapped.update(row)
    out.append(mapped)
  return out

def map_variant_to_MaveDB_scores(matches, mapped_info, row, extra):
  """
  Map variant information to a MaveDB score entry.

  Extracts the HGVS expression from mapped_info (using the first expression value).
  If no Variant Recoder matches are provided (matches is None), it calls join_information
  to directly join the variant details. Otherwise, it uses match_information to match the HGVS.
  """
  hgvs = mapped_info['expressions'][0]['value']

  if matches is None:
    return join_information(hgvs, mapped_info, row, extra)
  else:
    return match_information(hgvs, matches, row, extra)

def round_float_columns(row, round):
  """
  Round all float values in a row to a specified number of decimal places.
  """
  if round is not None:
    for i in row.keys():
      try:
        rounded = round(float(row[i]), round)
        row[i] = '{0:g}'.format(rounded)
      except:
        # If the conversion fails, leave the value unchanged
        pass
  return row

def map_scores_to_variants(scores, mappings, metadata, map_ids, matches=None, round=None, mapping_dict=None):
  """
  Map MaveDB scores to variant coordinates.

  For each score row:
    - Skip rows with special HGVS values (e.g. synonymous or wild-type) or missing scores.
    - Retrieve the corresponding mapping entry using the accession using mapping_dict.
    - Check for an URN mismatch between the score and mapping.
    - Round numeric values if requested.
    - If the mapping contains multiple members (phased variant), process each member;
      otherwise, process the mapping directly.

  Additional metadata (URN, publish_date, RefSeq, PubMed/DOI/URL) is merged into every output record.
  """

  # Extract reference sequence and publication identifiers
  refseq = None
  if len(metadata.get('targetGenes', [])) > 1:
    raise Exception("Multiple targets are not currently supported")
  else:
    for item in metadata.get('targetGenes', [{}])[0].get('externalIdentifiers', []):
      if item.get('identifier', {}).get('dbName') == 'RefSeq':
        refseq = item['identifier'].get('identifier')

  pubmed_list = []
  for pub in metadata.get('primaryPublicationIdentifiers', []):
    if pub.get('dbName') == 'PubMed':
      pubmed_list.append(pub.get('identifier'))
    else:
      warnings.warn("No PubMed ID found in metadata")
  pubmed = ",".join([p for p in pubmed_list if p])

  doi_list = []
  for pub in metadata.get('primaryPublicationIdentifiers', []):
    if pub.get('doi'):
      doi_list.append(pub.get('doi'))
    else:
      warnings.warn("No doi found in metadata")
  doi = ",".join([d for d in doi_list if d])

  url_list = []
  for pub in metadata.get('primaryPublicationIdentifiers', []):
    if pub.get('url'):
      url_list.append(pub.get('url'))
    else:
      warnings.warn("No URL found in metadata")
  url = ",".join([u for u in url_list if u])

  extra = {
    'urn'          : metadata.get('urn'),
    'publish_date' : metadata.get('experiment', {}).get('publishedDate'),
    'refseq'       : refseq,
    'pubmed'       : pubmed,
    'doi'          : doi,
    'url'          : url
  }

  out = []
  n_skipped_special = n_skipped_na = n_missing_map = n_urn_mismatch = n_no_post = 0

  for row in scores:
    subid = row.get('accession', 'na')

    # Skip rows with special HGVS values (e.g. synonymous, wild-type)
    if row.get('hgvs_pro') in ('_sy', '_wt', 'p.=') or row.get('hgvs_nt') in ('_sy', '_wt'):
      n_skipped_special += 1
      continue

    # Skip rows with missing score values
    if row.get('score') == "NA" or row.get('score') is None:
      n_skipped_na += 1
      continue

    # Retrieve the corresponding mapping entry using the pre-built mapping dictionary
    mapping = mapping_dict.get(subid)
    if mapping is None:
      log("accession_not_in_mappings", subid=subid)
      n_missing_map += 1
      continue

    # Check for URN mismatch between the score row and the mapping entry
    if 'mavedb_id' in mapping and subid != mapping['mavedb_id']:
      log("urn_mismatch_scores_vs_mappings", subid=subid, mapping_id=mapping.get('mavedb_id'))
      n_urn_mismatch += 1
      continue

    row = round_float_columns(row, round)

    # Some rows don't have post-mapped (i.e. no mapping, so only consider rows with mapping)
    if 'post_mapped' in mapping:
      mapped_info = mapping['post_mapped']
    else:
      log("no_post_mapped_in_mapping", subid=subid)
      n_no_post += 1
      continue

    # Process phased variants if multiple members exist; otherwise, process the mapping directly
    members = mapped_info.get("members", [])
    if members:
      for member in members:
        out += map_variant_to_MaveDB_scores(matches, member, row, extra)
    else:
      out += map_variant_to_MaveDB_scores(matches, mapped_info, row, extra)

  # Summary for this URN
  log("mapped_rows", n=len(out))
  log("skips_summary",
      skipped_special=n_skipped_special,
      skipped_missing_score=n_skipped_na,
      skipped_missing_mapping=n_missing_map,
      skipped_urn_mismatch=n_urn_mismatch,
      skipped_no_post_mapped=n_no_post)

  return out

def write_variant_mapping (f, map):
  """
  Write the final mapping between variants and MaveDB scores to an output TSV file.

  Constructs a header from the keys of the first output record (excluding unwanted fields),
  writes the header (with any 'hgvs_' prefixes removed), and then writes each record.
  """
  if map:
    with open(f, 'w') as csvfile:
      header = list(map[0].keys())
      header = [h for h in header if h not in ['HGVSp', 'index']]
      writer = csv.DictWriter(csvfile, delimiter="\t", fieldnames=header,
                              extrasaction='ignore')
      new_header = [h.replace('hgvs_', '') for h in header]
      header = OrderedDict(zip(header, new_header))
      writer.writerow(header)
      writer.writerows(map)
    log("write_done", out=f, rows=len(map))
    return True
  else:
    log("no_mappings_found", out=f)
    sys.exit(1)

if __name__ == "__main__":
  main()
