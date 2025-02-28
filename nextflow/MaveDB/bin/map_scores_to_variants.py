#!/usr/bin/env python3
import urllib.request
from urllib.parse import urlencode
import csv
import os
import json
import operator
import warnings
from collections import OrderedDict
from math import isclose
import re
import argparse

# Customise warning messages
def customshowwarning(message, category, filename, lineno, file=None, line=None):
    print("WARNING:", message)
warnings.showwarning = customshowwarning

def main():
  # Parse arguments
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

  # load MaveDB mappings, scores and HGVSP to variant matches
  print("Loading MaveDB data...", flush=True)
  scores = load_scores(args.scores)
  with open(args.mappings) as f:
    mappings = json.load(f)
  with open(args.metadata) as f:
    metadata = json.load(f)

  overhead = mappings[0]['id'] - 1
  mavedb_ids = [args.urn + "#" + str(i['id'] - overhead) for i in mappings ]
  print("DEBUG mavedb_ids: ", mavedb_ids)
  print("DEBUG overhead: ")
  
  if args.vr is not None:
    hgvsp2vars = load_vr_output(args.vr)
    print("DEBUG hgvsp2vars: ", hgvsp2vars)
  else:
    hgvsp2vars = None

  # create output file with variant location and respective scores
  print("Preparing mappings between variants and MaveDB scores...", flush=True)
  map = map_scores_to_variants(scores, mappings, metadata, mavedb_ids, hgvsp2vars, args.round)
  write_variant_mapping(args.output, map)

  print("Done: MaveDB score mapped to variants!", flush=True)
  return True

def load_vr_output (f):
  """Load Variant Recoder output"""
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
  """Load HGVSp to variant matches"""
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
  """Load MaveDB scores"""
  scores = []
  with open(f) as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
      scores.append(row)
  return scores

chrom = None
def get_chromosome (hgvs):
  """Lookup chromosome name in Ensembl REST API (unless already stored)"""
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

def join_information (hgvs, mapped_info, row, extra):
  """Join variant and MaveDB score information for a given HGVS"""
  
  print("DEBUG join_information: ")
  
  var = mapped_info['variation']
  ref = mapped_info['vrs_ref_allele_seq']
  
  print("var: ", var)
  print("ref: ", ref)

  mapped = OrderedDict(
    [("chr",   get_chromosome(hgvs)),
     ("start", var["location"]["interval"]["start"]["value"] + 1),
     ("end",   var["location"]["interval"]["end"]["value"]),
     ("ref",   ref),
     ("alt",   var["state"]["sequence"]),
     ("hgvs",  hgvs)])
  mapped.update(extra)
  mapped.update(row)
  return [mapped]

def match_information (hgvs, matches, row, extra):
  """Match a given HGVS to join variant and MaveDB score information"""
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

def map_variant_to_MaveDB_scores (matches, mapped_info, row, extra):
  
  hgvs = mapped_info['expressions'][0]['value']
  
  print("DEBUG map_variant_to_MaveDB_scores")
  print("len(mapped_info['expressions'])=", len(mapped_info['expressions']))
  print("mapped_info['expressions']: ", mapped_info['expressions'])
  print("mapped_info['expressions'][0]", mapped_info['expressions'][0])
  
  if matches is None:
    # HGVS genomic coordinates
    return join_information(hgvs, mapped_info, row, extra)
  else:
    # HGVS protein matches
    return match_information(hgvs, matches, row, extra)

def round_float_columns(row, round):
  """Round all values that are float"""
  if round is not None:
    for i in row.keys():
      try:
        rounded = round(float(row[i]), round)
        # Avoid rounding integers
        row[i] = '{0:g}'.format(rounded)
      except:
        # Not a number
        pass
  return row

def map_scores_to_variants (scores, mappings, metadata, map_ids, matches=None, round=None):
  """Map MaveDB scores to variants"""

  refseq = None
  if len(metadata['targetGenes']) > 1:
    raise Exception("multiple targets are not currently supported")
  else:
    for item in metadata['targetGenes'][0]['externalIdentifiers']:
      if item['identifier']['dbName'] == 'RefSeq':
        refseq = item['identifier']['identifier']

  pubmed_list = []
  for pub in metadata['primaryPublicationIdentifiers']:
    if pub['dbName'] == 'PubMed':
      pubmed_list.append(pub['identifier'])
  pubmed = ",".join(pubmed_list)

  extra = {
    'urn'          : metadata['urn'],
    'publish_date' : metadata['experiment']['publishedDate'],
    'refseq'       : refseq,
    'pubmed'       : pubmed,
  }

  out = []
  for row in scores:
    # Skip special HGVSp
    if row['hgvs_pro'] in ('_sy', '_wt', 'p.=') or row['hgvs_nt'] in ('_sy', '_wt'):
      continue

    # Skip missing values
    if row['score'] == "NA" or row['score'] is None:
      continue

    # Map available information
    if row['accession'] in map_ids:
      mapping = mappings[ map_ids.index(row['accession']) ]
    else:
      warnings.warn(row['accession'] + " not in mappings file")
      continue

    # check if accession is expected between mapping and score files
    if 'mavedb_id' in mapping.keys() and row['accession'] != mapping['mavedb_id']:
       warnings.warn("URN mismatch: trying to match " + row['accession'] + " from scores file with " + mapping['mavedb_id'] + " from mappings file")
       continue

    row = round_float_columns(row, round)
    mapped_info = mapping['postMapped']
    if mapped_info['type'] == "Haplotype":
      for member in mapped_info['members']:
        out += map_variant_to_MaveDB_scores(matches, member, row, extra)
    else:
      out += map_variant_to_MaveDB_scores(matches, mapped_info, row, extra)
  return out

def write_variant_mapping (f, map):
  """Write mapping between MaveDB scores and associated variants to output file"""

  with open(f, 'w') as csvfile:
    header = list(map[0].keys())
    header = [h for h in header if h not in ['HGVSp', 'index']]
    writer = csv.DictWriter(csvfile, delimiter="\t", fieldnames=header,
                            extrasaction='ignore')
    # prepare new header
    new_header = [h.replace('hgvs_', '') for h in header]
    header = OrderedDict(zip(header, new_header))

    writer.writerow(header)
    writer.writerows(map)
  return True

if __name__ == "__main__":
  main()