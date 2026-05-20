#!/usr/bin/env python3

"""
    Script to generate the gene annotation file for SpliceAI.
    This file contains data about the transcript that is being used for each gene.
    SpliceAI only annotates variants overlapping these transcripts.
    
    By passing --gff3, the supplied GFF3 is read in and transcript rows are retained. By default
    this keeps MANE_Select transcripts; with --gencode_primary it keeps gencode_primary
    protein-coding transcripts on main chromosomes.

    Each retained transcript is written as a separate SpliceAI annotation row. Exons are never
    merged across transcripts because that removes splice boundaries used by SpliceAI.

    Template provided by SpliceAI: https://github.com/Illumina/SpliceAI/blob/master/spliceai/annotations/grch38.txt

    Gene annotation file format:
        #NAME   CHROM   STRAND  TX_START    TX_END  EXON_START  EXON_END
        KRTAP27-1   21  -   30337013    30337694    30337013,   30337694,

    Options:
            --output_file   gene annotation output file (Optional. Default: gene_annotation.txt)
            --species       species name                (Optional. Default: homo_sapiens)
            --assembly      assembly version            (Optional. Default: 38)
            --gff3          GFF3 file path              (Mandatory)
            --gencode_primary   switch to filter for GENCODE primary instead of MANE Select
            --name_format       gene, transcript, or gene_transcript for the SpliceAI NAME field
"""

import argparse
import sys
import gzip

def fetch_transcripts_gff3(gff3_path, use_gencode_primary, name_format):
    annotations = {}
    open_func = gzip.open if gff3_path.endswith(".gz") else open
    tag_to_keep = "gencode_primary" if use_gencode_primary else "mane_select"
    transcript_features = {"transcript", "mRNA"}
    main_chromosomes = {str(chrom) for chrom in range(1, 23)} | {"X", "Y"}

    def detect_transcript_features(path, tag):
        features = set()
        with open_func(path, "rt") as feature_handle:
            for line in feature_handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue
                feature = fields[2]
                attrs_lower = fields[8].lower()
                if tag and tag in attrs_lower:
                    features.add(feature)
        return features

    def strip_prefix(value):
        if not value:
            return None
        val = value.split(",")[0]
        return val.split(":", 1)[1] if ":" in val else val

    def clean_chrom(chrom):
        return chrom[3:] if chrom.startswith("chr") else chrom

    def annotation_name(gene_name, transcript_id):
        if name_format == "gene":
            return gene_name
        if name_format == "transcript":
            return transcript_id
        return f"{gene_name}:{transcript_id}"

    # For mane_select, make the allowed transcript feature list to whatever appears tagged in the GFF to capture all
    if not use_gencode_primary:
        detected_features = detect_transcript_features(gff3_path, tag_to_keep)
        if detected_features:
            transcript_features = detected_features

    with open_func(gff3_path, "rt") as handle:
        transcripts_keep = {}
        gene_meta = {}
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _source, feature, start, end, _score, strand, _phase, attrs = fields
            chrom = clean_chrom(chrom)
            if chrom not in main_chromosomes:
                continue

            attr_dict = {}
            for entry in attrs.split(";"):
                if "=" in entry:
                    k, v = entry.split("=", 1)
                    attr_dict[k] = v

            if feature == "gene":
                gene_id = strip_prefix(attr_dict.get("ID"))
                gene_name = attr_dict.get("Name") or attr_dict.get("gene_name") or gene_id
                if gene_id:
                    gene_meta[gene_id] = {"name": gene_name, "chr": chrom, "strand": strand}
                    for transcript_id, transcript_gene_id in transcripts_keep.items():
                        if transcript_gene_id != gene_id:
                            continue
                        annotations[transcript_id]["name"] = annotation_name(gene_name, transcript_id)
                        annotations[transcript_id]["chr"] = chrom
                        annotations[transcript_id]["strand"] = strand
                continue

            # In gencode_primary mode, ignore everything except transcript and exon rows.
            if use_gencode_primary and feature not in transcript_features and feature != "exon":
                continue

            if feature in transcript_features:
                attrs_lower = attrs.lower()
                if tag_to_keep and tag_to_keep not in attrs_lower:
                    continue
                biotype = attr_dict.get("biotype") or attr_dict.get("gene_biotype")
                if use_gencode_primary and biotype != "protein_coding":
                    continue
                transcript_id = strip_prefix(attr_dict.get("transcript_id") or attr_dict.get("ID"))
                gene_id = strip_prefix(attr_dict.get("Parent"))
                if not (transcript_id and gene_id):
                    continue
                transcripts_keep[transcript_id] = gene_id
                gene_info = gene_meta.get(gene_id, {})
                gene_name = gene_info.get("name", gene_id)
                annotations[transcript_id] = {
                    "name": annotation_name(gene_name, transcript_id),
                    "chr": gene_info.get("chr", chrom),
                    "strand": gene_info.get("strand", strand),
                    "start": int(start),
                    "end": int(end),
                    "exons": set()
                }
                continue

            if feature == "exon":
                parents_raw = attr_dict.get("Parent", "")
                for parent in parents_raw.split(","):
                    transcript_id = strip_prefix(parent)
                    if transcript_id not in transcripts_keep:
                        continue
                    annotations[transcript_id]["exons"].add((int(start), int(end)))

    # convert exon sets to sorted lists
    formatted = {}
    for transcript_id, data in annotations.items():
        exons_sorted = sorted(data["exons"], key=lambda p: (p[0], p[1]))
        if not exons_sorted:
            continue
        exons_start = [str(p[0]) for p in exons_sorted]
        exons_end = [str(p[1]) for p in exons_sorted]
        formatted[transcript_id] = {
            "name": data["name"],
            "chr": data["chr"],
            "strand": data["strand"],
            "start": data["start"],
            "end": data["end"],
            "exons_start": exons_start,
            "exons_end": exons_end
        }
    return formatted

def sanity_checks(transcripts_list):
    ok = {}
    fail_list = []
    reason_counts = {}
    warning_counts = {}

    for gene, data in transcripts_list.items():
        original_pairs = [(int(s), int(e)) for s, e in zip(data["exons_start"], data["exons_end"])]
        # sort exons by start to keep output ordered before validation
        pairs = sorted(original_pairs, key=lambda p: (p[0], p[1]))
        data["exons_start"] = [str(p[0]) for p in pairs]
        data["exons_end"] = [str(p[1]) for p in pairs]
        reasons = []
        warnings = []
        if not pairs:
            reasons.append("no_exons")

        check = 1

        # overall start > end
        if pairs and data["start"] >= data["end"]:
            check = 0
            reasons.append("span_start_ge_end")

        # exon start/end validity
        for exon_start, exon_end in zip(data["exons_start"], data["exons_end"]):
            if int(exon_start) > int(exon_end):
                check = 0
                reasons.append("exon_start_gt_end")
                break

        # detect overlaps after sorting (start <= previous end)
        prev_end = None
        for exon_start, exon_end in zip(data["exons_start"], data["exons_end"]):
            if prev_end is not None and int(exon_start) <= int(prev_end):
                check = 0
                reasons.append("overlap")
                break
            prev_end = exon_end

        if pairs and (int(data["start"]) > pairs[0][0] or int(data["end"]) < pairs[-1][1]):
            check = 0
            reasons.append("transcript_span_does_not_cover_exons")

        # detect original out-of-order (before sorting, start decreases)
        prev_start_orig = None
        for exon_start, exon_end in original_pairs:
            if prev_start_orig is not None and exon_start < prev_start_orig:
                warnings.append("out_of_order")
                break
            prev_start_orig = exon_start

        if check == 0 or reasons:
            fail_list.append((gene, ";".join(reasons) if reasons else "unknown"))
            for r in reasons:
                reason_counts[r] = reason_counts.get(r, 0) + 1
        else:
            ok[gene] = data
            for w in warnings:
                warning_counts[w] = warning_counts.get(w, 0) + 1

    total = len(transcripts_list)
    print(f"[spliceai_annotation_file] Sanity check: total={total} pass={len(ok)} fail={len(fail_list)}", file=sys.stderr)
    if reason_counts:
        parts = [f"{reason}={count}" for reason, count in sorted(reason_counts.items())]
        print(f"[spliceai_annotation_file] Fail reasons: {', '.join(parts)}", file=sys.stderr)
    if warning_counts:
        warning_labels = {
            "overlap": "exons overlap after sorting (check input ordering)",
            "out_of_order": "exons not in ascending order in source"
        }
        parts = []
        for reason, count in sorted(warning_counts.items()):
            label = warning_labels.get(reason, reason)
            parts.append(f"{reason}={count} [{label}]")
        print(f"[spliceai_annotation_file] Warnings: {', '.join(parts)}", file=sys.stderr)

    return ok, fail_list

def write_output(transcripts_list, output_file):
    # Write to output file
    with open(output_file, "w") as f:
        f.write("#NAME\tCHROM\tSTRAND\tTX_START\tTX_END\tEXON_START\tEXON_END\n")

        for gene, data in transcripts_list.items():
            name = data.get("name", gene)
            chr = data["chr"]
            strand = data["strand"]
            # SpliceAI annotation files store starts as 0-based and ends as 1-based
            start = int(data["start"]) - 1
            end = data["end"]
            exons_start = ",".join(str(int(exon_start) - 1) for exon_start in data["exons_start"])
            exons_end = ",".join(data["exons_end"])

            f.write(f"{name}\t{chr}\t{strand}\t{start}\t{end}\t{exons_start},\t{exons_end},\n")


def main():
    parser = argparse.ArgumentParser(description="Generate the gene annotation file for SpliceAI")
    parser.add_argument("-o", "--output_file",
                        default="gene_annotation.txt",
                        help="output file (default: gene_annotation.txt)")
    parser.add_argument("-sp", "--species",
                        default="homo_sapiens",
                        help="species (default: homo_sapiens)")
    parser.add_argument("-a", "--assembly",
                        default="38",
                        help="species assembly (default: 38)")
    parser.add_argument("-r", "--release", required=True)
    parser.add_argument("--gff3", required=True,
                        help="GFF3 file path (required)")
    parser.add_argument("--gencode_primary", action="store_true",
                        help="Filter GFF3 transcripts to tag=gencode_primary instead of MANE_Select")
    parser.add_argument("--name_format", choices=["gene", "transcript", "gene_transcript"],
                        help="Annotation NAME field format (default: gene_transcript with --gencode_primary, otherwise gene)")
    args = parser.parse_args()

    output_file = args.output_file
    species = args.species
    assembly = args.assembly
    release = args.release
    if species.lower() not in ["homo_sapiens", "human"]:
        parser.error("Only human is currently supported")
    name_format = args.name_format or ("gene_transcript" if args.gencode_primary else "gene")
    filter_label = "gencode_primary" if args.gencode_primary else "MANE_Select"
    print(f"[spliceai_annotation_file] file={args.gff3} | filter={filter_label} | name_format={name_format}", file=sys.stderr)
    transcripts_list = fetch_transcripts_gff3(args.gff3, args.gencode_primary, name_format)

    ok, fail = sanity_checks(transcripts_list)
    sorted_list = dict(sorted(ok.items(), key=lambda kv: kv[0]))
    write_output(sorted_list, output_file)
    if fail:
        fail_strings = [f"{g}({r})" for g, r in fail]
        print("Sanity checks failed for the following genes: ", (", ").join(fail_strings), file=sys.stderr)


if __name__ == '__main__':
    main()
