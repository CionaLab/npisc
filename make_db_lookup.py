# %%

import sys
import re
import csv
import argparse
from collections import defaultdict
import networkx as nx
from networkx.algorithms.traversal.breadth_first_search import bfs_tree

from ontology import (
    build_graph,
    node_to_leaves,
    get_node_name,
    group_leaves_by_time,
    get_linked_lists,
    parse_url,
)

# %%

g_territory, g_stage = build_graph("developmental_ontology.txt")

pattern_stages = re.compile(r"^Stage \d+[a-z]?")
STAGES = [
    (i, j, g_stage.nodes[i]["synonym"])
    for l in get_linked_lists(g_stage)
    for (i, j) in l
    if pattern_stages.match(j)
    if pattern_stages.match(j)
]

# %%
with open('stages.tsv', 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['Term', 'Stage', 'Name'])
    for s in STAGES:
        writer.writerow(s)

# %%

# Find the root node
root_node = [n for n, d in g_territory.out_degree() if d==0]

with open('territories.tsv', 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['Term', 'Name'])
    for r in root_node:
        t = bfs_tree(g_territory, r, reverse=True)

        # Print the tree breadth first
        for node in t:
            writer.writerow((node, g_territory.nodes[node]["name"],))

# %%
STAGES

# %%
list(g_territory.nodes.data("name"))
