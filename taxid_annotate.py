import sys

from Bio import SeqIO


def cross_reference(ref_records, component, alternates=[]):
    # Step through the ref_records and search for a match
    # Return the best match
    best_match = None
    # print(component)
    all_components = [component] + alternates
    for component in all_components:
        for record in ref_records:
            if component.lower() in record.description.lower().replace(",", ""):
                # print(component, record.description)
                if best_match is None:
                    best_match = record
                    # print("new best match", best_match.description)
                elif record.name.startswith("NC_"):
                    if best_match.name.startswith("NC_"):
                        if len(record.seq) > len(best_match.seq):
                            best_match = record
                            # print("new best match", best_match.description)
                    else:
                        best_match = record
                        # print("new best match", best_match.description)
                elif not best_match.name.startswith("NC_") and len(record.seq) > len(best_match.seq):
                    best_match = record
                    # print("new best match", best_match.description)

    return best_match


def main(reference, panel_components):
    # First parse the reference fasta fna file using biopython
    ref_records = list(SeqIO.parse(reference, "fasta"))

    # Next step through the panel components file and search for matches in the ref_records
    # If a match is found, write the record to a new file
    with open(panel_components, "r") as panel, open("target_annotations.tsv", "a") as output:
        for line in panel:
            line = line.strip().split("\t")
            component = line[0]
            try:
                alternates = line[1].split(",")
            except IndexError:
                alternates = []

            # if line.startswith("#"):
            #     continue

            best_match = cross_reference(ref_records, component, alternates)

            if best_match is not None:
                ambig = str(best_match.seq).count("N")
                output.write(
                    f"{best_match.id}\t{component}\t{best_match.description}\t{len(best_match.seq)}\t{ambig}\n"
                )
            else:
                # for alternate in alternates:
                #     best_match = cross_reference(ref_records, alternate)
                #     if best_match is not None:
                #         ambig = str(best_match.seq).count("N")
                #         output.write(
                #             f"{best_match.id}\t{component}\t{best_match.description}\t{len(best_match.seq)}\t{ambig}\n"
                #         )
                #         break
                # if best_match is None:
                output.write("No match found\t" + component + "\n")


if __name__ == "__main__":
    # Get arguments from command line if they're present
    if len(sys.argv) == 1:
        reference = "genomic.fna"
        query = "panel_components.txt"
    elif len(sys.argv) == 2:
        reference = sys.argv[1]
        panel_components = "panel_components.txt"
    elif len(sys.argv) == 3:
        reference = sys.argv[1]
        panel_components = sys.argv[2]
    else:
        print("Error: Too many arguments!")
    main(reference, panel_components)
