#!/usr/bin/env python3

"""
    Script to generate the gene annotation file for SpliceAI.
    This file contains data about the transcript that is being used for each gene.
    SpliceAI only annotates variants overlapping these transcripts.

    Template provided by SpliceAI: https://github.com/Illumina/SpliceAI/blob/master/spliceai/annotations/grch38.txt

    Gene annotation file format:
        #NAME   CHROM   STRAND  TX_START    TX_END  EXON_START  EXON_END
        KRTAP27-1   21  -   30337013    30337694    30337013,   30337694,

    Options:
            --output_file   gene annotation output file (Optional. Default: gene_annotation.txt)
            --species       species name                (Optional. Default: homo_sapiens)
            --assembly      assembly version            (Optional. Default: 38)
            --host          core database host          (Mandatory)
            --port          host port                   (Mandatory)
            --user          host user                   (Mandatory)
            --release       core database version       (Mandatory)
"""

import argparse
import mysql.connector
from mysql.connector import Error


def fetch_transcripts(species, assembly, release, host, port, user):
    gene_annotation = {}
    database = f"{species}_core_{release}_{assembly}"

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
                    WHERE t.stable_id like 'ENST%' and t.biotype = 'protein_coding' and ga.attrib_type_id = 4
                """
    # For human select 'MANE_Select' transcripts
    if species == "homo_sapiens":
        sql_select = f"{sql_select} and atr.code = 'MANE_Select'"
    else:
        sql_select = f"{sql_select} and atr.code = 'is_canonical'"

    # Sort results
    sql_select = f"{sql_select} order by ga.value,s.name,t.seq_region_start,t.seq_region_end,e.seq_region_start,e.seq_region_end"

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

def sanity_checks(transcripts_list):
    fail_list = []

    for gene, data in transcripts_list.items():
        check = 1

        # transcript start > end
        if data["start"] >= data["end"]:
            check = 0
        
        # the exon start has to be bigger than the last exon end
        i = 0
        for exon_start in data["exons_start"]:
            if i > 0:
                last_exon_end = data["exons_end"][i - 1]
                if(int(exon_start) < int(last_exon_end)):
                    check = 0

            i += 1

        if check == 0:
            fail_list.append(gene)

    return fail_list

def write_output(transcripts_list, output_file):
    # Write to output file
    with open(output_file, "w") as f:
        f.write("#NAME\tCHROM\tSTRAND\tTX_START\tTX_END\tEXON_START\tEXON_END\n")

        for gene, data in transcripts_list.items():
            chr = data["chr"]
            strand = data["strand"]
            start = data["start"]
            end = data["end"]
            exons_start = (",").join(data["exons_start"])
            exons_end = (",").join(data["exons_end"])

            f.write(f"{gene}\t{chr}\t{strand}\t{start}\t{end}\t{exons_start}\t{exons_end}\n")


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
    parser.add_argument("--host", required=True)
    parser.add_argument("--port", required=True)
    parser.add_argument("--user", required=True)
    args = parser.parse_args()

    output_file = args.output_file
    species = args.species
    assembly = args.assembly
    release = args.release
    host = args.host
    port = args.port
    user = args.user

    transcripts_list = fetch_transcripts(species, assembly, release, host, port, user)

    check = sanity_checks(transcripts_list)

    if not check:
        write_output(transcripts_list, output_file)
    else:
        print("Sanity checks failed for the following genes: ", (", ").join(check))


if __name__ == '__main__':
    main()