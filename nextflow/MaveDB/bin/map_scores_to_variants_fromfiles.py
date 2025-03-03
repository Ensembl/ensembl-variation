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

# Customize warning messages to simply print the warning message
def customshowwarning(message, category, filename, lineno, file=None, line=None):
    print("WARNING:", message)
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
  print("Loading MaveDB data...", flush=True)
  
  scores = load_scores(args.scores)
  
  # Check if the scores file is empty, if so, print an error and exit - don't process this URN
  # This will occur because the previous step filtered out bad URN IDs i.e. "tmp:"
  if not scores:
    print(f"ERROR: The scores file '{args.scores}' for URN '{args.urn}' is empty (probably due to invalid URN IDs). Exiting.")
    sys.exit(1)
  
  with open(args.mappings) as f:
    mappings = json.load(f)
    
  with open(args.metadata) as f:
    metadata = json.load(f)
  
  # Extract MaveDB IDs from each mapping entry
  mavedb_ids = [i["mavedb_id"] for i in mappings["mapped_scores"]]
  
  # Pre-build a dictionary mapping accession IDs to mapping records
  # This makes lookups robust and order-independent
  mapping_dict = {mapping["mavedb_id"]: mapping for mapping in mappings["mapped_scores"]}
  
  # If a Variant Recoder output file is provided, load it; otherwise, set matches to None
  if args.vr is not None:
    hgvsp2vars = load_vr_output(args.vr)
  else:
    hgvsp2vars = None

  # Create the mapping between variant coordinates and MaveDB scores
  print("Preparing mappings between variants and MaveDB scores...", flush=True)
  mapped_data = map_scores_to_variants(scores, mappings, metadata, mavedb_ids, hgvsp2vars, args.round, mapping_dict)
  write_variant_mapping(args.output, mapped_data)

  print("Done: MaveDB score mapped to variants!", flush=True)
  return True

def load_vr_output (f):
  """
  Load Variant Recoder output.
  
  Iterates over each allele in the JSON data, skipping any entries that couldn't be parsed.
  For each valid allele, splits each VCF string (format: chr-start-ref-alt) and builds an
  ordered dictionary with variant details. Returns a dictionary mapping HGVS strings to lists of variants.
  """
  data = json.load(open(f))
  matches = {}
  for result in data:
    for allele in result:
      info = result[allele]
      if type(info) is list and "Unable to parse" in info[0]:
        continue
      hgvs = info["input"]
      for string in info["vcf_string"]:
        chr, start, ref, alt = string.split('-')
        end = int(start) + len(alt) - 1
        dict = OrderedDict([("HGVSp", hgvs),
                            ("chr",   chr),
                            ("start", start),
                            ("end",   end),
                            ("ref",   ref),
                            ("alt",   alt)])
        if hgvs not in matches:
          matches[hgvs] = []
        matches[hgvs].append(dict)
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
    for row in reader:
      scores.append(row)
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
  end = mapped_info["location"]["end"]
  # Extract the reference allele.
  ref = mapped_info["extensions"][0]["value"]
  # Extract the alternate allele.
  alt = mapped_info["state"]["sequence"]
  # Extract the HGVS expression.
  hgvs = mapped_info["expressions"][0]["value"]
            
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
  
  If the provided HGVS is not found in the matches, a warning is issued.
  For each matching entry, merge the match with extra metadata and the score row.
  """
  out = []
  if hgvs not in matches:
    warnings.warn(f"{hgvs} not found in HGVSp-variant matches")
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
        # Use the built-in round() function to round the value
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
  
  Additional metadata (URN, publish_date, RefSeq, PubMed) is merged into every output record.
  """
  
  refseq = None
  if len(metadata['targetGenes']) > 1:
    raise Exception("Multiple targets are not currently supported")
  else:
    for item in metadata['targetGenes'][0]['externalIdentifiers']:
      if item['identifier']['dbName'] == 'RefSeq':
        refseq = item['identifier']['identifier']

  pubmed_list = []
  for pub in metadata['primaryPublicationIdentifiers']:
    if pub['dbName'] == 'PubMed':
      pubmed_list.append(pub['identifier'])
    else:
      raise Exception("PubMed not found in metadata")
    
  pubmed = ",".join(pubmed_list)

  extra = {
    'urn'          : metadata['urn'],
    'publish_date' : metadata['experiment']['publishedDate'],
    'refseq'       : refseq,
    'pubmed'       : pubmed,
  }
  
  out = []
  for row in scores:
    # Skip rows with special HGVS values (e.g. synonymous, wild-type)
    if row['hgvs_pro'] in ('_sy', '_wt', 'p.=') or row['hgvs_nt'] in ('_sy', '_wt'):
      continue

    # Skip rows with missing score values
    if row['score'] == "NA" or row['score'] is None:
      continue

    # Retrieve the corresponding mapping entry using the pre-built mapping dictionary
    mapping = mapping_dict.get(row['accession'])
    if mapping is None:
        warnings.warn(row['accession'] + " not in mappings file")
        continue

    # Check for URN mismatch between the score row and the mapping entry
    if 'mavedb_id' in mapping.keys() and row['accession'] != mapping['mavedb_id']:
       warnings.warn("URN mismatch: trying to match " + row['accession'] + " from scores file with " + mapping['mavedb_id'] + " from mappings file")
       continue

    row = round_float_columns(row, round)
    mapped_info = mapping['post_mapped']
    
    # Process phased variants if multiple members exist; otherwise, process the mapping directly
    if mapped_info.get("members", []):
      for member in mapped_info['members']:
        out += map_variant_to_MaveDB_scores(matches, member, row, extra)
    else:
      out += map_variant_to_MaveDB_scores(matches, mapped_info, row, extra)
  
  return out

def write_variant_mapping (f, map):
  """
  Write the final mapping between variants and MaveDB scores to an output TSV file.
  
  Constructs a header from the keys of the first output record (excluding unwanted fields),
  writes the header (with any 'hgvs_' prefixes removed), and then writes each record.
  """
  with open(f, 'w') as csvfile:
    header = list(map[0].keys())
    header = [h for h in header if h not in ['HGVSp', 'index']]
    writer = csv.DictWriter(csvfile, delimiter="\t", fieldnames=header,
                            extrasaction='ignore')
    new_header = [h.replace('hgvs_', '') for h in header]
    header = OrderedDict(zip(header, new_header))
    writer.writerow(header)
    writer.writerows(map)
  return True

if __name__ == "__main__":
  main()