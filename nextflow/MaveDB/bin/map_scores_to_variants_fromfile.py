#!/usr/bin/env python3
import json
import csv
import argparse
import warnings
from collections import OrderedDict

# Custom warning messages
def customshowwarning(message, category, filename, lineno, file=None, line=None):
    print("WARNING:", message)
warnings.showwarning = customshowwarning

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description='Output file with MaveDB scores mapped to variants')
    parser.add_argument('--urn', type=str, help="MaveDB URN")
    parser.add_argument('--scores', type=str, help="Path to file with MaveDB URN scores")
    parser.add_argument('--mappings', type=str, help="Path to file with MaveDB URN mappings")
    parser.add_argument('--metadata', type=str, help="Path to file with MaveDB URN metadata")
    parser.add_argument('-o', '--output', type=str, help="Path to output file")
    parser.add_argument('--round', type=int, help="Number of decimal places for rounding values (default: not used)")
    args = parser.parse_args()

    # Load data
    print("Loading MaveDB data...", flush=True)
    scores = load_scores(args.scores)
    
    with open(args.mappings) as f:
        mappings = json.load(f)["mapped_scores"]  # Extract relevant section

    with open(args.metadata) as f:
        metadata = json.load(f)

    # Map MaveDB scores to variants
    print("Preparing mappings between variants and MaveDB scores...", flush=True)
    mapped_data = map_scores_to_variants(scores, mappings, metadata, args.urn, args.round)

    # Write output file
    write_variant_mapping(args.output, mapped_data)
    print("Done: MaveDB scores mapped to variants!", flush=True)

def load_scores(f):
    """Load MaveDB scores"""
    scores = []
    with open(f) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            scores.append(row)
    return scores

def extract_variant_info(mapping_entry):
    """Extract genomic variant information from mapping entry"""
    try:
        post_mapped = mapping_entry['post_mapped']
        members = post_mapped['members']

        variants = []
        for member in members:
            if member["type"] == "Allele":
                chrom = member['expressions'][0]['value'].split(':')[0]
                ref_allele = member['extensions'][0]['value']  # Extract reference allele
                alt_allele = member['state']['sequence']  # Extract alternative allele
                start = member['location']['start']
                end = member['location']['end']

                variants.append({
                    "chr": chrom,
                    "start": start + 1,  # Convert 0-based to 1-based coordinates
                    "end": end,
                    "ref": ref_allele,
                    "alt": alt_allele
                })

        return variants
    except KeyError as e:
        warnings.warn(f"Missing expected key in mapping entry: {e}")
        return None

def map_scores_to_variants(scores, mappings, metadata, urn, round=None):
    """Map MaveDB scores to genomic variants"""
    mapped_variants = []
    
    metadata_info = {
        "urn": metadata['urn'],
        "publish_date": metadata['experiment']['publishedDate'],
        "refseq": next((item['identifier']['identifier'] for item in metadata['targetGenes'][0]['externalIdentifiers'] if item['identifier']['dbName'] == 'RefSeq'), None),
        "pubmed": ",".join(pub['identifier'] for pub in metadata['primaryPublicationIdentifiers'] if pub['dbName'] == 'PubMed')
    }

    for row in scores:
        if row['hgvs_nt'] in ('_sy', '_wt') or row['score'] in ("NA", None):
            continue  # Skip synonymous and missing values

        matching_mapping = next((m for m in mappings if m['mavedb_id'] == row['accession']), None)
        if not matching_mapping:
            warnings.warn(f"Variant ID {row['accession']} not found in mappings file")
            continue

        variant_info_list = extract_variant_info(matching_mapping)
        if not variant_info_list:
            continue

        row = round_float_columns(row, round)
        row.update(metadata_info)

        for variant_info in variant_info_list:
            row_copy = row.copy()
            row_copy.update(variant_info)
            mapped_variants.append(row_copy)

    return mapped_variants

def round_float_columns(row, rounding):
    """Round all numerical values"""
    if round is not None:
        for key in row.keys():
            try:
                row[key] = '{0:g}'.format(round(float(row[key]), rounding))
            except ValueError:
                pass  # Not a number
    return row

def write_variant_mapping(output_file, mapped_data):
    """Write final mapped variants to output file"""
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ["chr", "start", "end", "ref", "alt", "urn", "publish_date", "refseq", "pubmed", "score"]
        writer = csv.DictWriter(csvfile, delimiter="\t", fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(mapped_data)

if __name__ == "__main__":
    main()