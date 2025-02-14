#!/usr/bin/env python3
import json
import argparse
import sys

def main(metadata_file, urn):
    # Load the JSON file
    with open(metadata_file, 'r') as f:
        try:
            data = json.load(f)
        except json.JSONDecodeError as e:
            print(f"Error decoding JSON: {e}")
            sys.exit(1)

    print(f"Extracting metadata for URN: {urn}")

    result = None
    # Check if the JSON is a dict with "experimentSets"
    if isinstance(data, dict) and "experimentSets" in data:
        # Loop over all experimentSets
        for exp_set in data.get("experimentSets", []):
            # Loop over all experiments in this experimentSet
            for experiment in exp_set.get("experiments", []):
                # Loop over each scoreSet in this experiment
                for score_set in experiment.get("scoreSets", []):
                    if score_set.get("urn") == urn:
                        result = score_set
                        break
                if result is not None:
                    break
            if result is not None:
                break
    # Fallback: if the JSON is a list of entries with top-level 'urn'
    elif isinstance(data, list):
        result = next((entry for entry in data if entry.get("urn") == urn), None)
    else:
        print("Unexpected JSON structure.")
        sys.exit(1)
    
    if result is None:
        print(f"URN '{urn}' not found in {metadata_file}")
        sys.exit(1)

    # Write the found metadata to a file called metadata.json
    with open("metadata.json", 'w') as f_out:
        json.dump(result, f_out, indent=2)

    print(f"Metadata for URN '{urn}' saved to metadata.json")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract metadata for a given URN from a JSON file"
    )
    parser.add_argument("--metadata_file", required=True,
                        help="Path to the large metadata JSON file")
    parser.add_argument("--urn", required=True,
                        help="Target URN to extract (e.g., 'urn:mavedb:00000001-a-1')")
    args = parser.parse_args()
    main(args.metadata_file, args.urn)

# TEST
# python extract_metadata.py --metadata_file /nfs/production/flicek/ensembl/variation/jma/maveDB-test/downloaded_data/main.json --urn "urn:mavedb:00000001-a-1"