import argparse
import datetime
import logging
import re
import sys

from Bio import SeqIO


def arg_parse():
    parser = argparse.ArgumentParser()
    input_group = parser.add_mutually_exclusive_group()
    input_group.add_argument("-f", "--fastq", help="Input fastq file, for single file filtering")
    input_group.add_argument(
        "-p", "--pair", nargs=2, help="Pair end reads in fastq format; R1 then R2. Final files will be synchronized."
    )
    parser.add_argument("-k", "--keep", nargs="*", help="Keep reads matching this regex", required=False)
    parser.add_argument("-d", "--drop", nargs="*", help="Drop reads matching this regex", required=False)
    parser.add_argument("-s", "--sync", help="Manual activation of syncronization", required=False)
    parser.add_argument(
        "-o",
        "--outfile",
        help="Specify a prefex name for your output files. Default is [Sample]_filtered.fastq.",
        required=False,
    )
    parser.add_argument(
        "-v",
        "--verbose",
        help="Creates logging file with information for debugging",
        required=False,
        action="store_true",
    )
    args = parser.parse_args()

    return args


def set_logging(args) -> None:
    if args.verbose:
        log_name = "fqfilter_{}.log".format(datetime.datetime.now().strftime("%y%m%d_%I_%M"))
        logging.basicConfig(filename=log_name, encoding="utf-8", level=getattr(logging, "DEBUG", None))
    else:
        logging.basicConfig(encoding="utf-8", level=getattr(logging, "INFO", None))
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(message)s")
    handler.setFormatter(formatter)
    root = logging.getLogger()
    root.addHandler(handler)
    logging.info("Logging started!")


def filter(file, keep, drop) -> list:
    total_input = 0
    dropped = 0
    pruned = []
    for record in SeqIO.parse(file, "fastq"):
        total_input += 1

        for reg in drop:
            if re.search(reg, str(record.seq)):
                dropped += 1
                break
        if len(keep):
            for reg in keep:
                if re.search(reg, str(record.seq)):
                    pruned.append(record)
                    break
        else:
            pruned.append(record)

    logging.info(f"Total Input Records: {total_input}")
    logging.info(f"Dropped Records: {dropped}")
    logging.info(f"Saved Records: {len(pruned)}")

    return pruned


def synchronize(read1, read2, file_prefix) -> None:
    r1_reads = {}
    for record in read1:
        read1_id = record.id.split(" ")
        r1_reads[read1_id[0]] = record

    r2_reads = {}
    for record in read2:
        read2_id = record.id.split(" ")
        r2_reads[read2_id[0]] = record

    r1_sync_file = open(f"{file_prefix}_R1.fastq", "w")
    r2_sync_file = open(f"{file_prefix}_R2.fastq", "w")

    for r2 in r2_reads.keys():
        if r1_reads.get(r2):
            SeqIO.write(r1_reads[r2], r1_sync_file, "fastq")
            SeqIO.write(r2_reads[r2], r2_sync_file, "fastq")


def main(args) -> None:
    keep_list = args.keep if args.keep else []
    drop_list = args.drop if args.drop else []

    if args.outfile:
        file_prefix = str(args.outfile)
    else:
        samplename = str(args.fastq).split(".", 1)[0]
        file_prefix = f"{samplename}_filtered"

    # - Debug: Print output file name
    logging.debug(f"Output name: {file_prefix}")
    # - Debug: Display input regex list
    logging.debug(f"Keep: {str(keep_list)}\nDrop: {str(drop_list)}")

    if args.pair:
        pruned_r1 = filter(args.pair[0], keep_list, drop_list)
        pruned_r2 = filter(args.pair[1], keep_list, drop_list)
        synchronize(pruned_r1, pruned_r2, file_prefix)

    else:
        pruned = filter(args.fastq, keep_list, drop_list)
        SeqIO.write(pruned, f"{file_prefix}.fastq", "fastq")


if __name__ == "__main__":
    args = arg_parse()
    set_logging(args)
    main(args)
