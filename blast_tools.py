import re
import argparse

'''
    blast_tools
    Chris Howard | 05/20/21
    
    Tools for manipulating/working with blast data
'''

## Remove Self Hits
# Given a blastn report file, save only the entries which are not self-hits
def self_ref_remove(in_lines):
    
    r_str = '^(chr.+:[0-9]+\-[0-9]+)\t(chr.+:[0-9]+\-[0-9]+)'
    out = []

    total = len(in_lines)
    prog = 0
    for line in in_lines:

        # Feedback
        current_prog = int((prog / total) * 100)
        if  current_prog % 5 == 0:
            print("Removing self hits.. " + str(current_prog) + "%", end="\r")
        prog += 1

        hit = re.match(r_str, line)
        if ( hit.group(1) == hit.group(2) ):
            next
        else: 
            out.append(str(line))
    
    print("Removing self hits.. Complete!")
    return out


def pick_best(in_lines):

    out = []

    return out


def main(args):

    data = []
    for line in open(args.input):
        data.append(line)

    if ( args.trimselfhits ):
        data = self_ref_remove(data)
    if ( args.pickbesthit ):
        data = pick_best(data)
    
    out = open(args.output, "w")
    for line in data:
        out.write(str(line))
    

if __name__ == "__main__":
    
    # Argument parsing
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-i", "--input", help="Input file", required=True)
    parser.add_argument("-v", "--verbose", help="Lots of status messages", action="store_true")
    parser.add_argument("-o", "--output", help="Designates a user-defined output file", default="output.txt")
    parser.add_argument("--trimselfhits", action="store_true")
    parser.add_argument("--pickbesthit", action="store_true")
    #parser.add_argument("-p", "--process", help="Selection of analysis tool", required=True)
    args = parser.parse_args()

    main(args)
