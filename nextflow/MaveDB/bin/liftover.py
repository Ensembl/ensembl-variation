#!/usr/bin/env python3
import os, json
from pyliftover import LiftOver

import argparse
parser = argparse.ArgumentParser(
                    description='Lift-over variants associated with MaveDB scores')
parser.add_argument('--mappings', type=str,
                    help="path to file with MaveDB mappings")
parser.add_argument('--mapped_variants', type=str,
                    help="path to file with variants mapped to MaveDB scores")
parser.add_argument('--reference', type=str, default="hg38",
                    help="genome to use as default (default: 'hg38')")
args = parser.parse_args()

reference = args.reference
mapped    = args.mapped_variants
mappings  = args.mappings

res = json.load(open(mappings))

def liftover_variants (mapped, genome, reference):
  # write file information with lifted-over coordinates to new file
  print(f"Converting coordinates from {genome} to {reference}...")
  chain = LiftOver(f"{genome}To{reference.capitalize()}.over.chain.gz")

  out = open(f"liftover_{mapped}", "w")
  with open(mapped) as f:
    header = f.readline()
    out.write(header)

    for line in f:
      l = line.split("\t")
      chr, start, end = l[0:3]

      conv = chain.convert_coordinate("chr" + chr, int(start))
      if len(conv) > 1:
        raise Exception("multiple coordinates returned")
      new_chr, new_start, strand, size = conv[0]

      conv = chain.convert_coordinate("chr" + chr, int(end))
      if len(conv) > 1:
        raise Exception("multiple coordinates returned")
      new_chr, new_end, strand, size = conv[0]

      l[0:3] = new_chr.replace("chr", ""), str(new_start), str(new_end)
      out.write('\t'.join(l))
  out.close()

if 'reference' in res['metadata']['extraMetadata']:
  genome = res['metadata']['extraMetadata']['reference']
else:
  genome = 'hg38'

if genome == reference:
  # just rename file if variants are already mapped to reference genome
  os.rename(mapped, f"liftover_{mapped}")
else:
  liftover_variants(mapped, genome, reference)
