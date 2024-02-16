#!/usr/bin/env python3
import argparse
import pandas as pd
import re
import os
import csv

parser = argparse.ArgumentParser()
parser.add_argument('--version', help="Release version", required=True)
parser.add_argument('--gtf', help="Assembly-specific GFF/GTF", required=True)
parser.add_argument('--outdir', default=".", help="Output directory")
parser.add_argument('--gene_symbols', default=None,
                    help="Lookup table with two columns (gene symbols and Ensembl identifiers) ")

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--go', help="GO terms annotation")
group.add_argument('--pheno', help="Phenotypes annotation")
args = parser.parse_args()

if args.go:
  plugin = 'GO'
  annot  = args.go
  ext    = 'gff'
  feat   = 'transcript'
elif args.pheno:
  plugin = 'phenotypes'
  annot  = args.pheno
  ext    = 'gvf'
  feat   = 'gene'

output = re.sub(r'(.*)-gca_(\d+)\.(\d+).*',
                f'\\1_gca\\2v\\3_{args.version}_VEP_{plugin}_plugin.{ext}',
                os.path.basename(args.gtf.lower()))

if not os.path.exists(args.outdir):
  os.makedirs(args.outdir)
output = args.outdir + "/" + output

# read assembly annotation
print(f"Preparing assembly annotation from {args.gtf}...", flush=True)
colnames = ['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
annot_pd = pd.read_csv(args.gtf, delimiter="\t", comment="#", header=None, names=colnames, dtype=str)
annot_pd = annot_pd[annot_pd['feature'].str.contains(feat)]
annot_pd = annot_pd.assign(gene=annot_pd['attribute'].str.extract(r'gene_id "(.*?)";'))
annot_pd = annot_pd.assign(gene_symbol=annot_pd['attribute'].str.extract(r'gene_name "(.*?)";'))
annot_pd = annot_pd.assign(transcript=annot_pd['attribute'].str.extract(r'transcript_id "(.*?)";'))

## read GO terms or Phenotypes annotation
print(f"Preparing {plugin} annotation from {annot}...", flush=True)
ref_pd = pd.read_csv(annot, delimiter="\t", comment="#", header=None, names=colnames, dtype=str)

if args.go:
  ref_pd = ref_pd.assign(gene_symbol=ref_pd['attribute'].str.extract(r'ID=(.*?);'))
elif args.pheno:
  ref_pd = ref_pd.assign(gene=ref_pd['attribute'].str.extract(r'id=(.*?);'))

# convert appropriate columns to numeric
annot_pd.start = annot_pd.start.astype(int)
annot_pd.end   = annot_pd.end.astype(int)
ref_pd.start   = ref_pd.start.astype(int)
ref_pd.end     = ref_pd.end.astype(int)

if args.gene_symbols is not None:
  # join with lookup table
  print(f"Preparing lookup table from {args.gene_symbols}...", flush=True)
  symbols = pd.read_csv(args.gene_symbols, delimiter="\t", names=['gene_symbol', 'gene'])
  ref_pd = pd.merge(ref_pd, symbols, on='gene', how='inner')

# join annotations based on gene symbols
print(f"Joining annotation...", flush=True)
joint = pd.merge(annot_pd, ref_pd, on='gene_symbol', how='inner')

if args.go:
  # prepare new assembly-specific GO annotations with backwards compatibiliy (i.e., transcript-based)
  joint = joint.assign(go_terms=joint['attribute_y'].str.extract(r';(Ontology_term=.*)'))
  joint = joint.assign(new_attribute="ID=" + joint['transcript'] + ';' + joint['go_terms'])
elif args.pheno:
  # prepare new assembly-specific Phenotypes annotations
  joint = joint.assign(phenotype=joint['attribute_y'].str.extract(r'; (phenotype=.*)'))
  joint = joint.assign(new_attribute="id=" + joint['gene_x'] + '; ' + joint['phenotype'])

# sort by genomic position
new_gtf = joint[['chr_x', 'source_y', 'feature_x', 'start_x', 'end_x',
                 'score_x', 'strand_x', 'frame_x', 'new_attribute']]
new_gtf = new_gtf.sort_values(by=['chr_x', 'start_x', 'end_x'])

if (len(new_gtf) == 0):
  raise Exception(f"ERROR: new pangenomes {plugin} annotation is empty (maybe no genes matched between annotations?)")  

# write to file
print(f"Writing new {plugin} annotation to {output}...", flush=True)
f = open(output, 'w')
if args.go:
  f.write('##gff-version 1.10\n')
new_gtf.to_csv(f, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)
f.close()
print(f"Done!", flush=True)
