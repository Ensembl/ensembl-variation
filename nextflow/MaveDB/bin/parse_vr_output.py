#!/usr/bin/env python3
import csv
import json

import argparse
parser = argparse.ArgumentParser(
                    description='Parse variant location from Variant Recoder')
parser.add_argument('filename', type=str,
                    help="path to Variant Recoder output")
parser.add_argument('-o', '--output', type=str, default="vr.txt",
                    help="path to output file")
args = parser.parse_args()

f = open(args.filename)
data = json.load(f)

rows = []
for result in data:
  for allele in result:
    info = result[allele]

    if type(info) is list and "Unable to parse" in info[0]:
      continue

    hgvsp  = info["input"]
    for string in info["vcf_string"]:
      chr, start, ref, alt = string.split('-')
      end = int(start) + len(alt) - 1
      rows.append( [hgvsp, chr, start, end, ref, alt] )
f.close()

w = open(args.output, 'w')
writer = csv.writer(w, delimiter='\t')
writer.writerow(["HGVSp", "chr", "start", "end", "ref", "alt"])
writer.writerows(rows)
w.close()
