#!/usr/bin/env python3

import argparse
import csv
import glob
import re
import sys
from pathlib import Path


MEMORY_TO_GB = {
    "KB": 1 / 1024 / 1024,
    "MB": 1 / 1024,
    "GB": 1,
    "TB": 1024,
}


def memory_to_gb(value):
    match = re.fullmatch(r"([0-9.]+)\s*([KMGT]B)", value.strip())
    if not match:
        raise ValueError("invalid memory value: {}".format(value))
    return float(match.group(1)) * MEMORY_TO_GB[match.group(2)]


def find_work_directory(work_dir, task_hash):
    prefix, suffix = task_hash.split("/", 1)
    matches = glob.glob(str(work_dir / prefix / (suffix + "*")))
    if len(matches) != 1:
        raise ValueError(
            "expected one work directory for {}, found {}".format(
                task_hash, len(matches)
            )
        )
    return Path(matches[0])


def parse_completed_tasks(trace_path):
    with trace_path.open() as trace_handle:
        next(trace_handle)
        for line_number, line in enumerate(trace_handle, start=2):
            fields = [field.strip() for field in line.rstrip("\n").split("\t")]
            try:
                process_index = fields.index("run_variant_recoder")
            except ValueError:
                continue

            # Locate fields relative to the process name so a recoverable joined
            # trace record does not discard a completed Variant Recoder task.
            if len(fields) <= process_index + 10:
                print(
                    "Skipping incomplete Variant Recoder trace row {}".format(
                        line_number
                    ),
                    file=sys.stderr,
                )
                continue

            if fields[process_index + 1] != "COMPLETED":
                continue

            yield {
                "hash": fields[process_index - 2],
                "peak_rss_gb": memory_to_gb(fields[process_index + 6]),
            }


def build_hints(trace_path, work_dir):
    hints = {}
    skipped = 0

    for task in parse_completed_tasks(trace_path):
        try:
            task_dir = find_work_directory(work_dir, task["hash"])
            command_log = (task_dir / ".command.log").read_text(errors="replace")
            command_script = (task_dir / ".command.sh").read_text(errors="replace")

            urn_match = re.search(r"MAVEDB_URN='([^']+)'", command_script)
            lines_match = re.search(r"hgvs_lines=(\d+)", command_log)
            if not urn_match or not lines_match:
                raise ValueError("URN or hgvs_lines missing from task files")

            urn = urn_match.group(1)
            hgvs_lines = int(lines_match.group(1))
            if hgvs_lines <= 0:
                raise ValueError("HGVS line count must be positive")

            previous = hints.get(urn)
            if previous is None or task["peak_rss_gb"] > previous["peak_rss_gb"]:
                hints[urn] = {
                    "hgvs_lines": hgvs_lines,
                    "peak_rss_gb": task["peak_rss_gb"],
                }
        except (OSError, ValueError) as error:
            skipped += 1
            print(
                "Skipping {}: {}".format(task["hash"], error),
                file=sys.stderr,
            )

    if not hints:
        raise RuntimeError("no completed Variant Recoder hints could be generated")

    return hints, skipped


def write_hints(hints, output_path):
    output_handle = sys.stdout if output_path is None else output_path.open("w")
    try:
        writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["urn", "hgvs_lines", "peak_rss_gb"])
        for urn in sorted(hints):
            hint = hints[urn]
            writer.writerow(
                [urn, hint["hgvs_lines"], "{:.3f}".format(hint["peak_rss_gb"])]
            )
    finally:
        if output_path is not None:
            output_handle.close()


def main():
    parser = argparse.ArgumentParser(
        description="Build per-URN Variant Recoder memory hints from a Nextflow trace"
    )
    parser.add_argument("--trace", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    hints, skipped = build_hints(args.trace, args.work_dir)
    write_hints(hints, args.output)
    print(
        "Wrote {} Variant Recoder memory hints; skipped {} tasks".format(
            len(hints), skipped
        ),
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
