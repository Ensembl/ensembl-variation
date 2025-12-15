#!/usr/bin/env python3

"""
    Script to generate the gene annotation file for SpliceAI.
    This file contains data about the transcript that is being used for each gene.
    SpliceAI only annotates variants overlapping these transcripts.

    Two routes:
      * DB mode (default): pulls MANE_Select protein-coding transcripts from the core DB (one per gene) and outputs their exon lists.
      * GFF3 mode (--gff3): reads the supplied GFF3, keeps protein-coding transcripts tagged MANE_Select (or gencode_primary when flagged), aggregates their exons per gene, and writes the same output format.

    Overlapping or nested exons are merged to the outermost span so the final exon list is non-overlapping.

    Template provided by SpliceAI: https://github.com/Illumina/SpliceAI/blob/master/spliceai/annotations/grch38.txt

    Gene annotation file format:
        #NAME   CHROM   STRAND  TX_START    TX_END  EXON_START  EXON_END
        KRTAP27-1   21  -   30337013    30337694    30337013,   30337694,

    Options:
            --output_file   gene annotation output file (Optional. Default: gene_annotation.txt)
            --species       species name                (Optional. Default: homo_sapiens)
            --assembly      assembly version            (Optional. Default: 38)
            --host          core database host          (Mandatory unless using --gff3)
            --port          host port                   (Mandatory unless using --gff3)
            --user          host user                   (Mandatory unless using --gff3)
            --database      override DB name            (Optional. Default: <species>_core_<release>_<assembly>)
            --release       core database version       (Mandatory)
            --gff3          GFF3 file path              (Optional. Enables file mode)
            --gencode_primary   switch to filter for GENCODE primary instead of MANE Select (GFF3 mode only)
"""

import argparse
import mysql.connector
from mysql.connector import Error
import sys
import gzip


def fetch_transcripts(species, assembly, release, host, port, user, database_name=None):
    gene_annotation = {}
    database = database_name or f"{species}_core_{release}_{assembly}"

    # For human we select 'MANE_Select' transcripts and the gene name is a gene attrib
    if species == "homo_sapiens":
        sql_select = """
                        SELECT ga.value,s.name,t.seq_region_strand,t.seq_region_start,t.seq_region_end,
                        e.seq_region_start,e.seq_region_end FROM transcript t
                        JOIN transcript_attrib ta ON t.transcript_id = ta.transcript_id
                        JOIN attrib_type atr ON ta.attrib_type_id = atr.attrib_type_id
                        JOIN seq_region s ON t.seq_region_id = s.seq_region_id
                        JOIN gene g ON g.gene_id = t.gene_id
                        JOIN gene_attrib ga ON g.gene_id = ga.gene_id
                        JOIN exon_transcript et ON t.transcript_id = et.transcript_id
                        JOIN exon e ON e.exon_id = et.exon_id
                        WHERE t.stable_id like 'ENST%' and t.biotype = 'protein_coding'
                        and ga.attrib_type_id = 4 and atr.code = 'MANE_Select'
                        order by ga.value,s.name,t.seq_region_start,t.seq_region_end,e.seq_region_start,e.seq_region_end
                """
    else:
        # For other species we select the canonical transcripts and the gene name is in xref
        sql_select = """
                        SELECT DISTINCT g.stable_id,s.name,t.seq_region_strand,t.seq_region_start,t.seq_region_end,
                        e.seq_region_start,e.seq_region_end FROM transcript t
                        JOIN transcript_attrib ta ON t.transcript_id = ta.transcript_id
                        JOIN attrib_type atr ON ta.attrib_type_id = atr.attrib_type_id
                        JOIN seq_region s ON t.seq_region_id = s.seq_region_id
                        JOIN gene g ON g.gene_id = t.gene_id
                        JOIN exon_transcript et ON t.transcript_id = et.transcript_id
                        JOIN exon e ON e.exon_id = et.exon_id
                        JOIN xref xr ON g.display_xref_id = xr.xref_id
                        WHERE t.stable_id like 'ENS%' and t.biotype = 'protein_coding' and atr.code = 'is_canonical'
                        order by xr.display_label,s.name,t.seq_region_start,t.seq_region_end,e.seq_region_start,e.seq_region_end
                     """

    connection = mysql.connector.connect(host=host,
                                         database=database,
                                         user=user,
                                         password='',
                                         port=port)

    try:
        if connection.is_connected():
            cursor = connection.cursor()
            cursor.execute(sql_select)
            data = cursor.fetchall()
            for row in data:
                strand = row[2]
                if strand == 1:
                    strand = "+"
                else:
                    strand = "-"

                if row[0] not in gene_annotation:
                    exons_start = []
                    exons_end = []
                    exons_start.append(str(row[5]))
                    exons_end.append(str(row[6]))

                    gene_annotation[row[0]] = {
                        "chr": row[1],
                        "strand": strand,
                        "start": row[3],
                        "end": row[4],
                        "exons_start": exons_start,
                        "exons_end": exons_end
                    }
                else:
                    gene_annotation[row[0]]["exons_start"].append(str(row[5]))
                    gene_annotation[row[0]]["exons_end"].append(str(row[6]))

    except Error as e:
        print(f"Error while connecting to MySQL {database}", e)
    finally:
        if connection.is_connected():
            cursor.close()
            connection.close()

    return gene_annotation

def fetch_transcripts_gff3(gff3_path, use_gencode_primary):
    gene_annotation = {}
    open_func = gzip.open if gff3_path.endswith(".gz") else open
    tag_to_keep = "gencode_primary" if use_gencode_primary else "mane_select"

    def strip_prefix(value):
        if not value:
            return None
        val = value.split(",")[0]
        return val.split(":", 1)[1] if ":" in val else val

    with open_func(gff3_path, "rt") as handle:
        transcripts_keep = {}
        gene_meta = {}
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, source, feature, start, end, score, strand, phase, attrs = fields
            chrom = chrom.replace("chr", "")

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
                    if gene_id in gene_annotation:
                        gene_annotation[gene_id]["name"] = gene_name
                        gene_annotation[gene_id]["chr"] = chrom
                        gene_annotation[gene_id]["strand"] = strand
                continue

            if feature not in ["transcript", "mRNA"] and feature != "exon":
                continue

            if feature in ["transcript", "mRNA"]:
                attrs_lower = attrs.lower()
                if tag_to_keep not in attrs_lower:
                    continue
                biotype = attr_dict.get("biotype") or attr_dict.get("gene_biotype")
                if biotype != "protein_coding":
                    continue
                transcript_id = strip_prefix(attr_dict.get("transcript_id") or attr_dict.get("ID"))
                gene_id = strip_prefix(attr_dict.get("Parent"))
                if not (transcript_id and gene_id):
                    continue
                transcripts_keep[transcript_id] = gene_id
                if gene_id not in gene_annotation:
                    gene_info = gene_meta.get(gene_id, {})
                    gene_annotation[gene_id] = {
                        "name": gene_info.get("name", gene_id),
                        "chr": gene_info.get("chr", chrom),
                        "strand": gene_info.get("strand", strand),
                        "exons": set()
                    }
                continue

            if feature == "exon":
                parents_raw = attr_dict.get("Parent", "")
                for parent in parents_raw.split(","):
                    transcript_id = strip_prefix(parent)
                    gene_id = transcripts_keep.get(transcript_id)
                    if not gene_id:
                        continue
                    if gene_id not in gene_annotation:
                        gene_info = gene_meta.get(gene_id, {})
                        gene_annotation[gene_id] = {
                            "name": gene_info.get("name", gene_id),
                            "chr": gene_info.get("chr", chrom),
                            "strand": gene_info.get("strand", strand),
                            "exons": set()
                        }
                    gene_annotation[gene_id]["exons"].add((int(start), int(end)))

    # convert exon sets to sorted lists and set start/end spans
    formatted = {}
    for gene_id, data in gene_annotation.items():
        exons_sorted = sorted(data["exons"], key=lambda p: (p[0], p[1]))
        if not exons_sorted:
            continue
        exons_start = [str(p[0]) for p in exons_sorted]
        exons_end = [str(p[1]) for p in exons_sorted]
        formatted[gene_id] = {
            "name": data["name"],
            "chr": data["chr"],
            "strand": data["strand"],
            "start": exons_sorted[0][0],
            "end": exons_sorted[-1][1],
            "exons_start": exons_start,
            "exons_end": exons_end
        }
    return formatted

def merge_overlapping_exons(pairs):
    """
    Collapse overlapping exon pairs, keeping the outermost span.
    Input pairs are tuples of ints sorted by start.
    Returns merged pairs (as strings) and a flag indicating if any merges happened.
    """
    merged = []
    merged_flag = False
    for start, end in pairs:
        if not merged:
            merged.append([start, end])
            continue
        prev_start, prev_end = merged[-1]
        if start <= prev_end:
            merged_flag = True
            if end > prev_end:
                merged[-1][1] = end
            # if the new exon is contained within the previous one, drop it
            continue
        merged.append([start, end])
    return [(str(p[0]), str(p[1])) for p in merged], merged_flag

def sanity_checks(transcripts_list):
    ok = {}
    fail_list = []
    reason_counts = {}
    warning_counts = {}
    warning_genes = {}

    for gene, data in transcripts_list.items():
        original_pairs = [(int(s), int(e)) for s, e in zip(data["exons_start"], data["exons_end"])]
        # sort exons by start to keep output ordered before validation
        pairs = sorted(original_pairs, key=lambda p: (p[0], p[1]))
        merged_pairs, had_merges = merge_overlapping_exons(pairs)
        data["exons_start"] = [p[0] for p in merged_pairs]
        data["exons_end"] = [p[1] for p in merged_pairs]
        reasons = []
        warnings = []
        if not pairs:
            reasons.append("no_exons")
        else:
            # span from first to last exon
            data["start"] = int(data["exons_start"][0])
            data["end"] = int(data["exons_end"][-1])

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
                warnings.append("overlap")
                break
            prev_end = exon_end

        if had_merges:
            warnings.append("overlap_merged")

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
                warning_genes.setdefault(w, []).append(gene)

    total = len(transcripts_list)
    print(f"[spliceai_annotation_file] Sanity check: total={total} pass={len(ok)} fail={len(fail_list)}", file=sys.stderr)
    if reason_counts:
        parts = [f"{reason}={count}" for reason, count in sorted(reason_counts.items())]
        print(f"[spliceai_annotation_file] Fail reasons: {', '.join(parts)}", file=sys.stderr)
    if warning_counts:
        warning_labels = {
            "overlap": "exons overlap after sorting (check input ordering)",
            "overlap_merged": "overlapping/nested exons merged to outer span",
            "out_of_order": "exons not in ascending order in source"
        }
        parts = []
        for reason, count in sorted(warning_counts.items()):
            label = warning_labels.get(reason, reason)
            parts.append(f"{reason}={count} [{label}]")
        print(f"[spliceai_annotation_file] Warnings: {', '.join(parts)}", file=sys.stderr)
        for reason, genes in warning_genes.items():
            sample = ", ".join(genes[:10])
            extra = warning_labels.get(reason, reason)

    return ok, fail_list

def write_output(transcripts_list, output_file):
    # Write to output file
    with open(output_file, "w") as f:
        f.write("#NAME\tCHROM\tSTRAND\tTX_START\tTX_END\tEXON_START\tEXON_END\n")

        for gene, data in transcripts_list.items():
            name = data.get("name", gene)
            chr = data["chr"]
            strand = data["strand"]
            start = data["start"]
            end = data["end"]
            exons_start = (",").join(data["exons_start"])
            exons_end = (",").join(data["exons_end"])

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
    parser.add_argument("--host")
    parser.add_argument("--port")
    parser.add_argument("--user")
    parser.add_argument("--database",
                        help="Override DB name (default: <species>_core_<release>_<assembly>)")
    parser.add_argument("--gff3",
                        help="GFF3 file path (enables file mode)")
    parser.add_argument("--gencode_primary", action="store_true",
                        help="Filter GFF3 transcripts to tag=gencode_primary instead of MANE_Select")
    args = parser.parse_args()

    output_file = args.output_file
    species = args.species
    assembly = args.assembly
    release = args.release
    host = args.host
    port = args.port
    user = args.user

    if args.gff3:
        if species.lower() not in ["homo_sapiens", "human"]:
            parser.error("GFF3 mode currently supports human only")
        filter_label = "gencode_primary" if args.gencode_primary else "MANE_Select"
        print(f"[spliceai_annotation_file] Mode: GFF3 | file={args.gff3} | filter={filter_label}", file=sys.stderr)
        transcripts_list = fetch_transcripts_gff3(args.gff3, args.gencode_primary)
    else:
        if not (host and port and user):
            parser.error("DB mode requires --host, --port and --user")
        if args.gencode_primary:
            print("[spliceai_annotation_file] Warning: --gencode_primary is ignored in DB mode (only MANE_Select is used)", file=sys.stderr)
        print(f"[spliceai_annotation_file] Mode: DB | db={args.database or f'{species}_core_{release}_{assembly}'} | host={host} | port={port} | user={user}", file=sys.stderr)
        transcripts_list = fetch_transcripts(species, assembly, release, host, port, user, args.database)

    ok, fail = sanity_checks(transcripts_list)
    sorted_list = dict(sorted(ok.items(), key=lambda kv: kv[0]))
    write_output(sorted_list, output_file)
    if fail:
        fail_strings = [f"{g}({r})" for g, r in fail]
        print("Sanity checks failed for the following genes: ", (", ").join(fail_strings), file=sys.stderr)


if __name__ == '__main__':
    main()
