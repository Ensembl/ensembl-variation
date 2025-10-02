#!/usr/bin/env python3
import json
import argparse
import sys
import os
import datetime

URN  = os.environ.get("MAVEDB_URN", "na")
STEP = os.environ.get("STEP", "extract_metadata")
def log(reason, subid="na", **kv):
    ts = datetime.datetime.now().astimezone().isoformat()
    extras = " ".join(f"{k}={json.dumps(v, ensure_ascii=False)}" for k, v in kv.items())
    sys.stderr.write(f"[{ts}][MaveDB][URN={URN}][STEP={STEP}][REASON={reason}][SUBID={subid}] {extras}\n")
    sys.stderr.flush()

def main(metadata_file, urn):
    # Load the JSON file
    try:
        with open(metadata_file, 'r') as f:
            data = json.load(f)
    except json.JSONDecodeError as e:
        log("json_decode_error", err=str(e), file=metadata_file)
        sys.exit(1)
    except FileNotFoundError:
        log("metadata_file_missing", file=metadata_file)
        sys.exit(1)

    log("extract_start", target_urn=urn)

    selected_entry = None
    for experiment_set in data.get("experimentSets", []):
        for experiment in experiment_set.get("experiments", []):
            for score_set in experiment.get("scoreSets", []):
                if score_set.get("urn") == urn:
                    selected_entry = score_set
                    break

    # Reformat the extracted data to match the downstream-expecting JSON
    if selected_entry:
        log("metadata_found", title=selected_entry.get("title", ""), shortName=selected_entry.get("license", {}).get("shortName", ""))
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
        log("metadata_not_found", file=metadata_file, target_urn=urn)
        sys.stderr.flush()
        sys.exit(1)

    # Save the formatted data
    with open("metadata.json", "w") as outfile:
        json.dump(formatted_data, outfile, indent=4)
    log("metadata_written", out="metadata.json", bytes=os.path.getsize("metadata.json"))

    # Write LICENCE.txt
    short = formatted_data.get('license', {}).get('shortName', '')
    with open("LICENCE.txt", "w") as f:
        f.write(short)
    try:
        sz = os.path.getsize("LICENCE.txt")
    except OSError:
        sz = 0
    log("licence_written", out="LICENCE.txt", shortName=short, bytes=sz)

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