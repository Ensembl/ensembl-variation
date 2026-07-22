#!/usr/bin/env python3
"""Inspect HGVS expressions in MaveDB public-dump mapping files."""

import argparse
import json
import sys

HGVS_PRIORITY = ("hgvs.p", "hgvs.g")


def load_records(path):
    with open(path) as handle:
        records = json.load(handle)
    if not isinstance(records, list):
        raise ValueError("Mappings JSON must be a top-level array of mapped variant records")
    return records


def iter_expressions(vrs_obj):
    if not isinstance(vrs_obj, dict):
        return

    for expression in vrs_obj.get("expressions") or []:
        if isinstance(expression, dict):
            yield expression

    for member in vrs_obj.get("members") or []:
        yield from iter_expressions(member)


def iter_current_postmapped(records):
    for record in records:
        if not isinstance(record, dict):
            continue
        if record.get("current") is False:
            continue
        post_mapped = record.get("postMapped")
        if post_mapped:
            yield post_mapped


def collect_values(records, syntax=None):
    values = []
    for post_mapped in iter_current_postmapped(records):
        for expression in iter_expressions(post_mapped):
            expr_syntax = expression.get("syntax")
            value = expression.get("value")
            if value and (syntax is None or expr_syntax == syntax):
                values.append(value)
    return sorted(set(values))


def command_type(args):
    records = load_records(args.mappings)
    syntaxes = {
        expression.get("syntax")
        for post_mapped in iter_current_postmapped(records)
        for expression in iter_expressions(post_mapped)
        if expression.get("syntax") in HGVS_PRIORITY
    }

    for syntax in HGVS_PRIORITY:
        if syntax in syntaxes:
            print(syntax)
            return 0

    print("No current postMapped HGVS expression found", file=sys.stderr)
    return 2


def command_extract(args):
    records = load_records(args.mappings)
    for value in collect_values(records, args.syntax):
        print(value)
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    type_parser = subparsers.add_parser("type", help="print the mapping HGVS type used for pipeline branching")
    type_parser.add_argument("mappings")
    type_parser.set_defaults(func=command_type)

    extract_parser = subparsers.add_parser("extract", help="print sorted unique HGVS expression values")
    extract_parser.add_argument("--syntax", required=True, choices=HGVS_PRIORITY)
    extract_parser.add_argument("mappings")
    extract_parser.set_defaults(func=command_extract)

    args = parser.parse_args()
    try:
        return args.func(args)
    except BrokenPipeError:
        try:
            sys.stdout.close()
        except OSError:
            pass
        return 0
    except Exception as exc:
        print(str(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
