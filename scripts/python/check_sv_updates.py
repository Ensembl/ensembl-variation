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

def get_studies_db(release, var_host, var_port, var_user):

    studies_from_db = {}

    database_name = f"homo_sapiens_variation_{release}_38"
    connection = mysql.connector.connect(host=var_host,
                                         database=database_name,
                                         user=var_user,
                                         password='',
                                         port=var_port)

    sql_query_select = "select st.name from study st left join source so on so.source_id = st.source_id where so.name = 'dbVar'"

    try:
        if connection.is_connected():
            cursor = connection.cursor()
            cursor.execute(sql_query_select)
            data = cursor.fetchall()
            print ("Fetching studies from database...")
            for row in data:
                studies_from_db[row[0]] = 1
            print ("Fetching studies from database... done")

    except Error as e:
        print("Error while connecting to MySQL", e)
    finally:
        if connection.is_connected():
            cursor.close()
            connection.close()
            print("MySQL connection is closed")

    return studies_from_db


def main():
    parser = argparse.ArgumentParser(description="Checks structural variant studies "
                                                 "to be updated/imported ")
    parser.add_argument("-sp", "--species",
                        default="Homo_sapiens",
                        help="species (default: Homo_sapiens)")
    parser.add_argument("-a", "--assembly",
                        default="GRCh38",
                        help="species assembly (default: GRCh38)")
    parser.add_argument("-f", "--format",
                        default="gvf")
    parser.add_argument("-r", "--release")
    parser.add_argument("--host")
    parser.add_argument("--port")
    parser.add_argument("--user")
    args = parser.parse_args()

    species = args.species
    assembly = args.assembly
    format = args.format
    release = args.release
    host = args.host
    port = args.port
    user = args.user
    files_dir = f"/pub/dbVar/data/{species}/by_study/{format}"

    ftp = FTP(DBVAR_HOST)
    ftp.login()
    ftp.cwd(files_dir)

    # get list of studies from database
    current_studies = get_studies_db(release, host, port, user)

    # use the date from production db
    # first update the import script to write to production db
    today = date.today()
    current_year = today.year

    out = []
    ftp.retrlines('LIST', out.append)

    # store info for each study id
    file_list = {}

    # output files
    dest_dir = os.getcwd()
    output_file_update = os.path.join(dest_dir, "studies_to_update.txt")
    output_file_notdb = os.path.join(dest_dir, "studies_to_import.txt")

    # specific to dbVar - create separate function
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

        if (assembly == file_assembly):
            file_list[file_study] = str(file_date)

    with open(output_file_update, "w") as f:
        f.write("Study name\tDate last updated on dbVar FTP\n")
        for st in file_list.keys():
            if st in current_studies:
                f.write(st + "\t" + file_list[st] + "\n")

if __name__ == '__main__':
    main()
