#!/usr/bin/env python3
import urllib.request
from urllib.parse import urlencode
import csv, os, json, operator, warnings
from collections import OrderedDict
from math import isclose
import re, argparse, sys, datetime

warnings.resetwarnings()
warnings.filterwarnings("always")

URN  = os.environ.get("MAVEDB_URN", "na")
STEP = os.environ.get("STEP", "map_scores")
def log(reason, subid="na", **kv):
    ts = datetime.datetime.now().astimezone().isoformat()
    extra = " ".join(f"{k}={json.dumps(v, ensure_ascii=False)}" for k,v in kv.items())
    sys.stderr.write(f"[{ts}][MaveDB][URN={URN}][STEP={STEP}][REASON={reason}][SUBID={subid}] {extra}\n")
    sys.stderr.flush()

# Pass Python warnings through uniform logger
def customshowwarning(message, category, filename, lineno, file=None, line=None):
    log("warning", msg=str(message))
warnings.showwarning = customshowwarning

def main():
  parser = argparse.ArgumentParser(description='Output file with MaveDB scores mapped to variants')
  parser.add_argument('--vr', type=str, help="path to file containg Variant Recoder output with 'vcf_string' enabled (optional)")
  parser.add_argument('--urn', type=str, help="MaveDB URN")
  parser.add_argument('--scores', type=str, help="path to file with MaveDB URN scores")
  parser.add_argument('--mappings', type=str, help="path to file with MaveDB URN mappings")
  parser.add_argument('--metadata', type=str, help="path to file with MaveDB URN metadata")
  parser.add_argument('-o', '--output', type=str, help="path to output file")
  parser.add_argument('--round', type=int, help="Number of decimal places for rounding values (default: not used)")
  args = parser.parse_args()

  # Load the scores, mappings, metadata
  log("load_start", scores=args.scores, mappings=args.mappings, metadata=args.metadata, vr=args.vr)

  scores = load_scores(args.scores)
  if not scores:
    log("scores_empty", target=args.scores, urn=args.urn)
    sys.exit(1)
  if len(scores) < 50:
    log("scores_small", n=len(scores), urn=args.urn)

  with open(args.mappings) as f:
    mappings = json.load(f)
  with open(args.metadata) as f:
    metadata = json.load(f)

  # Pre-build dict
  mapped_scores = mappings.get("mapped_scores", [])
  mavedb_ids = [i.get("mavedb_id") for i in mapped_scores if "mavedb_id" in i]
  mapping_dict = {m["mavedb_id"]: m for m in mapped_scores if "mavedb_id" in m}
  log("mappings_loaded", n=len(mapped_scores), unique_ids=len(mapping_dict))

  # Optional VR
  if args.vr is not None:
    hgvsp2vars = load_vr_output(args.vr)
  else:
    hgvsp2vars = None

  log("map_prepare")
  mapped_data = map_scores_to_variants(scores, mappings, metadata, mavedb_ids, hgvsp2vars, args.round, mapping_dict)
  write_variant_mapping(args.output, mapped_data)

  log("done", out=args.output)
  return True

def load_vr_output(f):
  data = json.load(open(f))
  # If any warnings, surface them; if all warnings => exit
  if any("warnings" in item for item in data):
    log("vr_has_warnings", file=f)
    if all("warnings" in item for item in data):
      log("vr_only_unable_to_parse", file=f)
      sys.exit(1)

  matches = {}
  for result in data:
    for allele in result:
      info = result[allele]
      if isinstance(info, list) and ("Unable to parse" in info[0] or "skipped" in info[0]):
        # mark the specific HGVSp that failed as SUBID (allele key is the input)
        log("vr_unable_to_parse", subid=allele, msg=info[0] if info else "")
        continue
      hgvs = info["input"]
      for string in info["vcf_string"]:
        chr, start, ref, alt = string.split('-')
        end = int(start) + len(alt) - 1
        d = OrderedDict([("HGVSp", hgvs), ("chr", chr), ("start", start), ("end", end), ("ref", ref), ("alt", alt)])
        matches.setdefault(hgvs, []).append(d)
  log("vr_loaded", n=sum(len(v) for v in matches.values()))
  return matches

def load_HGVSp_to_variant_matches(f):
  matches = {}
  with open(f) as csvfile:
    reader = csv.DictReader(csvfile, delimiter="\t")
    for row in reader:
      hgvs = row['HGVSp']
      matches.setdefault(hgvs, []).append(row)
  return matches

def load_scores(f):
  scores = []
  with open(f) as csvfile:
    reader = csv.DictReader(csvfile)
    reader.fieldnames = [field.strip() for field in reader.fieldnames]
    for row in reader:
      clean_row = {k: (v.strip() if isinstance(v, str) else v) for k, v in row.items()}
      scores.append(clean_row)
  log("scores_loaded", n=len(scores))
  return scores

chrom = None
def get_chromosome(hgvs):
  global chrom
  if (chrom is None):
    chrom = hgvs.split(":")[0]
    url = f"https://rest.ensembl.org/info/assembly/homo_sapiens/{chrom}?"
    data = urlencode({"synonyms": 1, "content-type": "application/json"})
    res = urllib.request.urlopen(url + data).read()
    res = json.loads(res)
    chrom = [each["name"] for each in res["synonyms"] if each['dbname'] == "UCSC"][0]
    chrom = chrom.replace("chr", "")
  return chrom

def join_information(hgvs, mapped_info, row, extra):
  start = mapped_info["location"]["start"]
  end   = mapped_info["location"]["end"]
  ref   = mapped_info["extensions"][0]["value"]
  alt   = mapped_info["state"]["sequence"]
  hgvs  = mapped_info["expressions"][0]["value"]

  mapped = OrderedDict([
    ("chr", get_chromosome(hgvs)),
    ("start", start + 1),
    ("end", end),
    ("ref", ref),
    ("alt", alt),
    ("hgvs", hgvs)
  ])
  mapped.update(extra); mapped.update(row)
  return [mapped]

def match_information(hgvs, matches, row, extra):
  out = []
  if hgvs not in matches:
    # per-HGVSp miss
    log("hgvsp_not_in_vr_matches", subid=hgvs)
    return out
  for match in matches[hgvs]:
    mapped = match
    mapped['hgvs'] = match['HGVSp']
    mapped.update(extra); mapped.update(row)
    out.append(mapped)
  return out

def map_variant_to_MaveDB_scores(matches, mapped_info, row, extra):
  hgvs = mapped_info['expressions'][0]['value']
  if matches is None:
    return join_information(hgvs, mapped_info, row, extra)
  else:
    return match_information(hgvs, matches, row, extra)

def round_float_columns(row, round):
  if round is not None:
    for i in row.keys():
      try:
        rounded = round(float(row[i]), round)
        row[i] = '{0:g}'.format(rounded)
      except:
        pass
  return row

def map_scores_to_variants(scores, mappings, metadata, map_ids, matches=None, round=None, mapping_dict=None):
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
      warnings.warn("No PubMed ID found in metadata")
  pubmed = ",".join(pubmed_list)

  doi_list = []
  for pub in metadata['primaryPublicationIdentifiers']:
    if pub['doi']:
      doi_list.append(pub['doi'])
    else:
      warnings.warn("No doi found in metadata")
  doi = ",".join(doi_list)

  url_list = []
  for pub in metadata['primaryPublicationIdentifiers']:
    if pub['url']:
      url_list.append(pub['url'])
    else:
      warnings.warn("No URL found in metadata")
  url = ",".join(url_list)

  extra = {
    'urn'          : metadata['urn'],
    'publish_date' : metadata['experiment']['publishedDate'],
    'refseq'       : refseq,
    'pubmed'       : pubmed,
    'doi'          : doi,
    'url'          : url
  }

  out = []
  n_skipped_special = n_skipped_na = n_missing_map = n_urn_mismatch = n_no_post = 0

  for row in scores:
    subid = row.get('accession', 'na')

    # special HGVS rows
    if row.get('hgvs_pro') in ('_sy', '_wt', 'p.=') or row.get('hgvs_nt') in ('_sy', '_wt'):
      n_skipped_special += 1
      continue

    # missing score
    if row.get('score') == "NA" or row.get('score') is None:
      n_skipped_na += 1
      continue

    # lookup mapping by accession
    mapping = mapping_dict.get(subid)
    if mapping is None:
        log("accession_not_in_mappings", subid=subid)
        n_missing_map += 1
        continue

    # accession mismatch vs mapping
    if 'mavedb_id' in mapping and subid != mapping['mavedb_id']:
       log("urn_mismatch_scores_vs_mappings", subid=subid, mapping_id=mapping['mavedb_id'])
       n_urn_mismatch += 1
       continue

    row = round_float_columns(row, round)

    if 'post_mapped' in mapping:
      mapped_info = mapping['post_mapped']
    else:
      log("no_post_mapped_in_mapping", subid=subid)
      n_no_post += 1
      continue

    # phased vs single
    members = mapped_info.get("members", [])
    if members:
      for member in members:
        out += map_variant_to_MaveDB_scores(matches, member, row, extra)
    else:
      out += map_variant_to_MaveDB_scores(matches, mapped_info, row, extra)

  # summary for this URN
  log("mapped_rows", n=len(out))
  log("skips_summary",
      skipped_special=n_skipped_special,
      skipped_missing_score=n_skipped_na,
      skipped_missing_mapping=n_missing_map,
      skipped_urn_mismatch=n_urn_mismatch,
      skipped_no_post_mapped=n_no_post)

  return out

def write_variant_mapping(f, map):
  if map:
    with open(f, 'w') as csvfile:
      header = list(map[0].keys())
      header = [h for h in header if h not in ['HGVSp', 'index']]
      writer = csv.DictWriter(csvfile, delimiter="\t", fieldnames=header, extrasaction='ignore')
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