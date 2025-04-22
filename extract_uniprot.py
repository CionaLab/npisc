import argparse
import sys

import pandas as pd

# Latest release at:
# https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.dat.gz

# State machine constants
STATE_INITIAL = 0
STATE_AC = 1
STATE_DE = 2
STATE_DR = 3

if __name__ == "__main__":
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

    with args.input as input_file:
        data = []
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
                    _, go_acc, go_text, go_source = line.split("; ")
                    go.append(f"{go_acc} {go_text} {go_source[:-1]}")
                elif line.startswith("//"):
                    state = STATE_INITIAL
                    data.append([accession_number, full_name, go])

        df = pd.DataFrame(data, columns=["uniprot", "fullname", "go"])
        df["go"] = df["go"].str.join("; ")
        df.to_csv(args.output, index=False)
