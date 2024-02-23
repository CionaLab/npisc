import sys
import re
import csv
import argparse
from collections import defaultdict

from ontology import (
    build_graph,
    node_to_leaves,
    get_node_name,
    group_leaves_by_time,
    get_linked_lists,
    parse_url,
)


def main():
    parser = argparse.ArgumentParser(description="Process ontology and tsv file.")
    parser.add_argument(
        "--obo",
        "-b",
        type=str,
        required=True,
        help="Path to the ontology OBO file",
    )
    parser.add_argument(
        "--input",
        "-i",
        type=argparse.FileType("r"),
        default=sys.stdin,
        help="Path to the input TSV file",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="Path to the output TSV file",
    )
    args = parser.parse_args()

    g_territory, g_stage = build_graph(args.obo)

    pattern_stages = re.compile(r"^Stage \d+[a-z]?")
    STAGES = [
        (i, j)
        for l in get_linked_lists(g_stage)
        for (i, j) in l
        if pattern_stages.match(j)
        if pattern_stages.match(j)
    ]
    DICT_STAGES = {k: i for i, (k, _) in enumerate(STAGES)}
    DICT_ID_STAGES = {v: i for i, (_, v) in enumerate(STAGES)}

    map_territory = group_leaves_by_time(g_territory, node_to_leaves(g_territory))

    pattern_blast = re.compile(r"^[AaBb]\d+\.\d+\*?$")

    d3 = defaultdict(lambda: defaultdict(list))
    for (id, start, end), v in map_territory.items():
        for i, _ in enumerate(STAGES):
            if i <= DICT_STAGES.get(end, 0) and i >= DICT_STAGES.get(
                start, len(STAGES)
            ):
                d3[i][get_node_name(g_territory, id)].extend(
                    [
                        name
                        for n in v
                        if pattern_blast.match(name := get_node_name(g_territory, n))
                    ]
                )
    visited = set()

    with args.input as fr, args.output as fw:
        reader = csv.DictReader(fr, delimiter="\t")
        writer = csv.DictWriter(fw, fieldnames=reader.fieldnames + ["Territory_eq"])
        writer.writeheader()
        for row in reader:
            params = parse_url(row["URL"])
            territory = row["Territory"]
            gene = row["Gene"]
            stage = row["Stage"]
            if (params["biomaterial_id"], stage, territory, gene) not in visited:
                visited.add((params["biomaterial_id"], stage, territory, gene))
                stage = pattern_stages.match(row["Stage"]).group(0)
                for t_eq in d3[DICT_ID_STAGES[stage]][row["Territory"]]:
                    writer.writerow({**row, "Territory_eq": t_eq})


if __name__ == "__main__":
    main()
