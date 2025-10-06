#!/usr/bin/env python3

import argparse
import json
import re
import sys
import urllib.parse
import urllib.request
import os
from ftplib import FTP, error_perm
import pandas as pd
import io


"""
Script to download files from the OpenTargets FTP and prepare them for the
phenotype pipeline.

Usage:
    python download_cgc_file.py
"""


HOST = "ftp.ebi.ac.uk"
BASE_DIR = "/pub/databases/opentargets/platform"
EVIDENCE_DIR = "output/association_by_datasource_direct"

def find_parquet_files(ftp, pathname):
    current = ftp.pwd()
    try:
        ftp.cwd(pathname)
    except error_perm as e:
        print ("Error: " + str(e))
        return  # file or private directory

    for name in ftp.nlst():
        if name.endswith(".parquet"):
            yield f"{pathname}/{name}"
        elif name == "_SUCCESS":
            continue
        else:
            yield from find_parquet_files(ftp, f"{pathname}/{name}")

    ftp.cwd(current)


def walk_ftp(host, dirname):
    ftp = FTP(host)
    ftp.login()

    try:
        for filename in find_parquet_files(ftp, dirname):
            url = f"ftp://{host}{filename}"

            with urllib.request.urlopen(url) as resp:
                data = resp.read()

            try:
                df = pd.read_parquet(io.BytesIO(data))
            except Exception as e:
                raise RuntimeError("Failed to read parquet file.")

            # Convert DataFrame to a list of JSON serialisable python dicts
            records = json.loads(df.to_json(orient='records', date_format='iso'))

            for rec in records:
                yield rec

    finally:
        ftp.quit()


def main():
    parser = argparse.ArgumentParser(description="Retrieve target - disease "
                                                 "evidences from the Open "
                                                 "Targets Platform")
    parser.add_argument("-r", "--release",
                        default="latest",
                        help="release version (default: latest)")
    parser.add_argument("-d", "--dest_dir",
                        default=os.getcwd())
    args = parser.parse_args()

    release = args.release
    out_dir = args.dest_dir
    evidence_dir = f"{BASE_DIR}/{release}/{EVIDENCE_DIR}"

    output_file = os.path.join(out_dir, f"cgc_input_{release}.json")

    if release != "latest":
        if not re.fullmatch(r"\d\d\.\d\d", release):
            parser.exit(status=2,
                        message=f"error: invalid release version: "
                                f"{args.release}\n")
        elif float(release) < 25.03:
            parser.exit(status=2,
                        message=f"error: release {args.release} "
                                f"not supported\n")

    sys.stderr.write(f"Fetching data from Open Target Platform "
                     f"(release: {release})\n")

    with open(output_file, "w") as f:
      for obj in walk_ftp(HOST, evidence_dir):
        try:
            disease_id = obj["diseaseId"]
        except KeyError:
            continue

        if obj["datasourceId"] != "cancer_gene_census":
            continue

        f.write(json.dumps(obj) + "\n")

if __name__ == '__main__':
    main()
