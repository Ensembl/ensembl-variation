"""
Version: 1.0 (2022-10-14)
"""

import argparse
import re
import sys
import os
from ftplib import FTP, error_perm
from datetime import date
from time import strptime
import mysql.connector
from mysql.connector import Error

DBVAR_HOST = "ftp.ncbi.nlm.nih.gov"

def get_studies_db(type, release, var_host, var_port, var_user, assembly, species):

    studies_from_db = {}

    # remove characters from the assembly - keep version
    # ex: GRCh38 -> 38
    db_assembly = re.sub('\D', '', assembly)

    species_name = species.lower()

    if type == "variation":
        database_name = f"{species_name}_variation_{release}_{db_assembly}"
        sql_query_select = """ SELECT st.name, st.study_id FROM study st 
                               LEFT JOIN source so ON so.source_id = st.source_id 
                               WHERE so.name = %s """
        params = ["dbVar"]
    else:
        database_name = type
        sql_query_select = """ SELECT cr.result, cr.date FROM check_result cr 
                              LEFT JOIN check_dictionary cd ON cd.check_id = cr.check_id 
                              LEFT JOIN ensvar_db en ON en.ensvardb_id = cr.ensvardb_id 
                              WHERE cd.description = 'structural_variation_study' AND cr.is_current = 1 
                              AND en.name LIKE %s """
        params = [f"{species_name}_variation_%_{db_assembly}"]

    connection = mysql.connector.connect(host=var_host,
                                         database=database_name,
                                         user=var_user,
                                         password='',
                                         port=var_port)

    try:
        if connection.is_connected():
            cursor = connection.cursor()
            cursor.execute(sql_query_select, params)
            data = cursor.fetchall()
            print (f"Fetching studies from {database_name}...")
            for row in data:
                studies_from_db[row[0]] = row[1]
            print (f"Fetching studies from {database_name}... done")

    except Error as e:
        print("Error while connecting to MySQL", e)
    finally:
        if connection.is_connected():
            cursor.close()
            connection.close()
            print("MySQL connection is closed")

    return studies_from_db


def main():
    parser = argparse.ArgumentParser(description="Structural Variant studies to be updated for a release")
    parser.add_argument("-sp", "--species",
                        default="Homo_sapiens",
                        help="species (default: Homo_sapiens)")
    parser.add_argument("-a", "--assembly",
                        default="GRCh38",
                        help="species assembly (default: GRCh38)")
    parser.add_argument("-f", "--format",
                        default="gvf")
    parser.add_argument("-r", "--release", required=True)
    parser.add_argument("--host", required=True)
    parser.add_argument("--prod-host", required=True)
    parser.add_argument("--port", required=True)
    parser.add_argument("--prod-port", required=True)
    parser.add_argument("--user", required=True)
    parser.add_argument("--print-all", action='store_true',
                        help=""" use option --print_all to return all studies from the variation db 
                        even the ones not found in production db """)
    args = parser.parse_args()

    species = args.species
    assembly = args.assembly
    format = args.format
    release = args.release
    host = args.host
    prod_host = args.prod_host
    port = args.port
    prod_port = args.prod_port
    user = args.user
    option_print = args.print_all
    files_dir = f"/pub/dbVar/data/{species}/by_study/{format}"

    ftp = FTP(DBVAR_HOST)
    ftp.login()
    ftp.cwd(files_dir)

    # get dictionary of studies from the Variation database
    current_studies = get_studies_db("variation", release, host, port, user, assembly, species)

    # get dictionary of studies from production db
    # this db contains the date last time the studies where imported into variation database
    studies_from_production = get_studies_db("production", release, prod_host, prod_port, user, assembly, species)

    # use the date from production db
    # first update the import script to write to production db
    current_year = date.today().year

    out = []
    ftp.retrlines('LIST', out.append)

    # store info for each study id
    file_list = {}

    # output file
    dest_dir = os.getcwd()
    output_file_update = os.path.join(dest_dir, "studies_to_update.txt")

    for line in out:
        line.strip()

        line_split = line.split()
        month = strptime(line_split[5],'%b').tm_mon
        day = line_split[6]
        more_det = line_split[7]
        if not re.fullmatch(r"\d\d\d\d", more_det):
            file_year = current_year
        else:
            file_year = more_det

        file_date = date(int(file_year), int(month), int(day))

        filename = line_split[8].replace(f".{format}.gz", "")
        filename_split = filename.split(".", 2);
        file_study = filename_split[0]
        file_assembly = filename_split[1]

        if assembly == file_assembly and file_study not in file_list:
            file_list[file_study] = file_date

    with open(output_file_update, "w") as f:
        f.write("Study name\tDate last updated on dbVar FTP\tComments\n")
        for st in file_list.keys():
            if option_print == False:
                if st in current_studies and st in studies_from_production and studies_from_production[st].date() < file_list[st]:
                    f.write(f"{st}\t{str(file_list[st])}\tPlease update study\n")
            else:
                if st in current_studies and st not in studies_from_production:
                    f.write(f"{st}\t{str(file_list[st])}\tStudy not found in production db\n")

if __name__ == '__main__':
    main()
