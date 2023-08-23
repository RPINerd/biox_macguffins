"""
    TSV/CSV to DB | RPINerd, 08/22/23

    Given either a tab- or comma-delimited file, build a simple SQLite database
"""

import argparse
import csv
import logging
import sqlite3


def main(args):
    # TODO handle multiple input files each as a table
    # Set up the database
    conn = sqlite3.connect(args.output)
    c = conn.cursor()

    # Create the table
    c.execute("DROP TABLE IF EXISTS data")
    c.execute(
        """CREATE TABLE data (
            accession TEXT PRIMARY KEY,
            taxid TEXT,
            gi TEXT
            )"""
    )

    # Read in the input file
    with open(args.input, "r") as input_file:
        if args.tab:
            reader = csv.reader(input_file, delimiter="\t")
        else:
            reader = csv.reader(input_file, delimiter=",")
        for row in reader:
            # Skip the header
            if row[0] == "accession":
                continue

            accession = row[1]
            taxid = row[2]
            gi = row[3]
            c.execute(
                "INSERT INTO data (accession, taxid, gi) VALUES (?, ?, ?)",
                (accession, taxid, gi),
            )

    # Commit the changes and close the connection
    conn.commit()
    conn.close()


if __name__ == "__main__":
    # Argument parsing
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-i", "--input", help="Input file", required=True)
    parser.add_argument("-o", "--output", help="Output database file", required=True)
    parser.add_argument("-v", "--verbose", help="Lots of status messages", action="store_true")
    parser.add_argument("-t", "--tab", help="Input file is tab-delimited", action="store_true")
    args = parser.parse_args()

    # Set up logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG, format="%(levelname)s: %(message)s")
    else:
        logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    main(args)
