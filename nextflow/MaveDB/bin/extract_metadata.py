#!/usr/bin/env python3
import json
import argparse
import sys
import os

def main(metadata_file, urn):
    # Load the JSON file
    with open(metadata_file, 'r') as f:
        try:
            data = json.load(f)
        except json.JSONDecodeError as e:
            print(f"Error decoding JSON: {e}")
            sys.exit(1)

    print(f"Extracting metadata for URN: {urn}")
    
    selected_entry = None
    for experiment_set in data.get("experimentSets", []):
        for experiment in experiment_set.get("experiments", []):
            for score_set in experiment.get("scoreSets", []):
                if score_set.get("urn") == urn:
                    selected_entry = score_set
                    break
                
    # Reformat the extracted data to match the json format expected later in the pipeline
    # This was a pragmatic approach so that the whole pipeline wasn't re-written
    # This is to cope with the fact that the pipeline was written for API -yielded json structures, 
    # which differ from data-dump download -yielded json structures
    if selected_entry:
        formatted_data = {
            "abstractText": selected_entry.get("abstractText", ""),
            "contributors": [],
            "createdBy": {
                "firstName": selected_entry["createdBy"].get("firstName", ""),
                "lastName": selected_entry["createdBy"].get("lastName", ""),
                "orcidId": selected_entry["createdBy"].get("orcidId", ""),
                "recordType": "User"
            },
            "creationDate": selected_entry.get("creationDate", ""),
            "datasetColumns": selected_entry.get("datasetColumns", {}),
            "doiIdentifiers": selected_entry.get("doiIdentifiers", []),
            "experiment": {
                "abstractText": selected_entry.get("abstractText", ""),
                "contributors": [],
                "createdBy": {
                    "firstName": selected_entry["createdBy"].get("firstName", ""),
                    "lastName": selected_entry["createdBy"].get("lastName", ""),
                    "orcidId": selected_entry["createdBy"].get("orcidId", ""),
                    "recordType": "User"
                },
                "creationDate": selected_entry.get("creationDate", ""),
                "doiIdentifiers": selected_entry.get("doiIdentifiers", []),
                "experimentSetUrn": experiment_set.get("urn"),
                "extraMetadata": selected_entry.get("extraMetadata", {}),
                "keywords": [],
                "methodText": selected_entry.get("methodText", ""),
                "modificationDate": selected_entry.get("modificationDate", ""),
                "modifiedBy": {
                    "firstName": selected_entry["modifiedBy"].get("firstName", ""),
                    "lastName": selected_entry["modifiedBy"].get("lastName", ""),
                    "orcidId": selected_entry["modifiedBy"].get("orcidId", ""),
                    "recordType": "User"
                },
                "primaryPublicationIdentifiers": selected_entry.get("primaryPublicationIdentifiers", []),
                "publishedDate": selected_entry.get("publishedDate", ""),
                "rawReadIdentifiers": selected_entry.get("rawReadIdentifiers", []),
                "recordType": "Experiment",
                "scoreSetUrns": [selected_entry.get("urn")],
                "secondaryPublicationIdentifiers": [],
                "shortDescription": selected_entry.get("shortDescription", ""),
                "title": selected_entry.get("title", ""),
                "urn": selected_entry.get("urn"),
            },
            "externalLinks": {},
            "extraMetadata": selected_entry.get("extraMetadata", {}),
            "license": {
                "active": True,
                "id": selected_entry.get("license", {}).get("id", 1),
                "link": selected_entry.get("license", {}).get("link", ""),
                "longName": selected_entry.get("license", {}).get("longName", ""),
                "recordType": "ShortLicense",
                "shortName": selected_entry.get("license", {}).get("shortName", ""),
                "version": selected_entry.get("license", {}).get("version", ""),
            },
            "mappingState": "complete",
            "metaAnalyzedByScoreSetUrns": [],
            "metaAnalyzesScoreSetUrns": [],
            "methodText": selected_entry.get("methodText", ""),
            "modificationDate": selected_entry.get("modificationDate", ""),
            "modifiedBy": {
                "firstName": selected_entry["modifiedBy"].get("firstName", ""),
                "lastName": selected_entry["modifiedBy"].get("lastName", ""),
                "orcidId": selected_entry["modifiedBy"].get("orcidId", ""),
                "recordType": "User"
            },
            "numVariants": selected_entry.get("numVariants", ""),
            "primaryPublicationIdentifiers": selected_entry.get("primaryPublicationIdentifiers", []),
            "private": selected_entry.get("private", ""),
            "processingState": selected_entry.get("processingState", ""),
            "publishedDate": selected_entry.get("publishedDate", ""),
            "recordType": "ScoreSet",
            "secondaryPublicationIdentifiers": [],
            "shortDescription": selected_entry.get("shortDescription", ""),
            "targetGenes": selected_entry.get("targetGenes", []),
            "title": selected_entry.get("title", ""),
            "urn": selected_entry.get("urn"),
    }
    else:
        print(f"ERROR: extract_metadata.py - no matching entry found for '{urn}' in metadata file '{metadata_file}'. Exiting.")
        sys.exit(1)

    # Save the formatted data
    with open("metadata.json", "w") as outfile:
        json.dump(formatted_data, outfile, indent=4)
        
    print(f"Metadata for URN '{urn}' saved to metadata.json")
    
    # Output a file containing the licence to allow downstream filtering based on this
    with open("LICENCE.txt", "w") as f:
        f.write(formatted_data['license']['shortName'])
    
    print(f"Licence for URN '{urn}' saved to LICENCE.txt")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract metadata for a given URN from a JSON file"
    )
    parser.add_argument("--metadata_file", required=True,
                        help="Path to the large metadata JSON file from the MaveDB data dump")
    parser.add_argument("--urn", required=True,
                        help="Target URN to extract (e.g., 'urn:mavedb:00000001-a-1')")
    args = parser.parse_args()
    main(args.metadata_file, args.urn)

## TEST
# python /hps/software/users/ensembl/variation/fairbrot/ensembl-variation/nextflow/MaveDB/bin/extract_metadata.py --metadata_file /nfs/production/flicek/ensembl/variation/jma/maveDB-test/mavedb_dbdump_data/main.json --urn "urn:mavedb:00001204-a-4"