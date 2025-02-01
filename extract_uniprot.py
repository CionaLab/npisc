import argparse
import csv
import sys

# State machine constants
STATE_INITIAL = 0
STATE_AC = 1
STATE_DE = 2
STATE_DR = 3

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Extract UniProt data")
parser.add_argument(
    "-i",
    "--input",
    default=sys.stdin,
    type=argparse.FileType("r"),
    help="Input file path (default: stdin)",
)
parser.add_argument(
    "-o",
    "--output",
    default=sys.stdout,
    type=argparse.FileType("w"),
    help="Output file path (default: stdout)",
)
args = parser.parse_args()

# Read UniProt records from the input file and write the data to the output file
with args.input as input_file, args.output as output_file:
    writer = csv.writer(output_file)
    writer.writerow(["uniprot", "fullname", "go"])

    state = STATE_INITIAL
    accession_number = ""
    full_name = ""
    go = []

    for line in input_file:
        if line.startswith("ID"):
            state = STATE_AC
            accession_number = ""
            full_name = ""
            go = []
        elif state == STATE_AC:
            if line.startswith("AC"):
                accession_number = line.split()[1].split(";")[0]
                state = STATE_DE
        elif state == STATE_DE:
            if line.startswith("DE   RecName:"):
                full_name = line.split("=")[1].split(";")[0]
                state = STATE_DR
        elif state == STATE_DR:
            if line.startswith("DR   GO; "):
                _, go_acc, go_text, _ = line.split("; ")
                go.append(f"{go_acc} {go_text}")
            elif line.startswith("//"):
                state = STATE_INITIAL
                writer.writerow([accession_number, full_name, "; ".join(go)])
