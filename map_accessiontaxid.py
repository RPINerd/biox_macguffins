"""
    Accession to TaxID Mapping | RPINerd, 08/22/23

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


def indexing(ref_path) -> dict:
    # Create a dictionary of the reference sequences
    ref_dict = {}
    priority_order = ["nucl_gb", "nucl_wgs", "nucl_wgs_extra", "dead_nucl", "dead_wgs"]
    for reference in priority_order:
        with open(ref_path + reference + ".a2t", "r") as ref_file:
            for record in ref_file:
                records = record.strip().split("\t")
                accession = records[1]
                taxid = records[2]
                gi = records[3]
                # If the accession number is already in the dictionary, skip it
                # This enforces the priority by saving the first unique accession number it finds
                if accession not in ref_dict:
                    ref_dict[accession] = [taxid, gi]

    return ref_dict


def accession_to_taxid(accession, ref_dict) -> str:
    # If the accession number is in the reference dictionary, return the taxid and gi
    if accession in ref_dict:
        return ref_dict[accession][0], ref_dict[accession][1]
    # If the accession number is not in the reference dictionary, return "NA"
    else:
        return "NA"


def main(args):
    # Create a dictionary of the reference sequences
    # ref_dict = indexing(args.reference)

    # Create a connection to the provided database
    conn = sqlite3.connect(args.reference)
    c = conn.cursor()

    # Create a list of the output sequences from stepping throught the input sequences and querying the database
    output_list = []
    for record in SeqIO.parse(args.input, "fasta"):
        accession = record.id
        # print(record.description)
        # print(f"Processing {accession}")
        # print(f"ID: {record.id}")
        info = record.description
        # print(f"Info: {info}")

        # taxid, gi = accession_to_taxid(record[0], ref_dict)
        # output_list.append([record[0], record[1], taxid, record[2]])
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
