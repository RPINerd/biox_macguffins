"""
    Accession to TaxID Mapping | RPINerd, 08/29/23

    Given a fasta input file and reference database, create a tsv output that maps the accession numbers to a taxid

    Priority for taxid lookup is:
      nucl_gb
      nucl_wgs
      nucl_wgs_extra
      dead_nucl
      dead_wgs
"""


import argparse
import sqlite3

from Bio import SeqIO


def main(args):
    # Create a connection to the provided database
    conn = sqlite3.connect(args.reference)
    c = conn.cursor()

    # Create a list of the output sequences from stepping throught the input sequences and querying the database
    output_list = []
    for record in SeqIO.parse(args.input, "fasta"):
        accession = record.id
        info = record.description

        c.execute("SELECT taxid, gi FROM data WHERE accession=?", (accession,))
        result = c.fetchone()
        if result is None:
            output_list.append([accession, "NA", info, record.seq])
        else:
            output_list.append([accession, result[0], info, str(record.seq)])

    # Write the output to a tsv file
    with open(args.output, "w") as output:
        # Header line
        output.write("accession\ttaxid\tinfo\tsequence\n")
        for i in range(len(output_list)):
            output.write("\t".join(output_list[i]) + "\n")


if __name__ == "__main__":
    # Parse the arguments
    parser = argparse.ArgumentParser(description="Map accession numbers to taxids")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input fasta file",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="Output tsv file",
    )
    parser.add_argument(
        "-r",
        "--reference",
        type=str,
        required=True,
        help="Reference database file",
    )
    args = parser.parse_args()

    main(args)
