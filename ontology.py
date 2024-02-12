from typing import List, Tuple, Optional, Dict
import re
import csv
import networkx as nx

STAGES = [
    ("CirobuD:0000021", "Stage 0"),
    ("CirobuD:0000022", "Stage 1"),
    ("CirobuD:0000023", "Stage 2"),
    ("CirobuD:0000024", "Stage 3"),
    ("CirobuD:0000025", "Stage 4"),
    ("CirobuD:0000026", "Stage 5a"),
    ("CirobuD:0000027", "Stage 5b"),
    ("CirobuD:0000028", "Stage 6a"),
    ("CirobuD:0000029", "Stage 6b"),
    ("CirobuD:0000030", "Stage 7"),
    ("CirobuD:0000031", "Stage 8"),
    ("CirobuD:0000032", "Stage 9"),
    ("CirobuD:0000033", "Stage 10"),
    ("CirobuD:0000034", "Stage 11"),
    ("CirobuD:0000035", "Stage 12"),
    ("CirobuD:0000036", "Stage 13"),
    ("CirobuD:0000037", "Stage 14"),
    ("CirobuD:0000038", "Stage 15"),
    ("CirobuD:0000039", "Stage 16"),
    ("CirobuD:0000040", "Stage 17"),
    ("CirobuD:0000041", "Stage 18"),
    ("CirobuD:0000042", "Stage 19"),
    ("CirobuD:0000043", "Stage 20"),
    ("CirobuD:0000044", "Stage 21"),
    ("CirobuD:0000045", "Stage 22"),
    ("CirobuD:0000046", "Stage 23"),
    ("CirobuD:0000047", "Stage 24"),
    ("CirobuD:0000048", "Stage 25"),
    ("CirobuD:0000049", "Stage 26"),
    ("CirobuD:0000050", "Stage 27"),
    ("CirobuD:0000051", "Stage 28"),
    ("CirobuD:0000052", "Stage 29"),
    ("CirobuD:0000053", "Stage 30"),
    ("CirobuD:0000054", "Stage 31"),
    ("CirobuD:0000055", "Stage 32"),
    ("CirobuD:0000056", "Stage 33"),
    ("CirobuD:0000057", "Stage 34"),
    ("CirobuD:0000058", "Stage 35"),
    ("CirobuD:0000059", "Stage 36"),
    ("CirobuD:0000060", "Stage 37"),
    ("CirobuD:0000061", "Stage 38"),
    ("CirobuD:0000062", "Stage 39"),
    ("CirobuD:0000063", "Stage 40"),
    ("CirobuD:0000064", "Stage 41"),
    ("CirobuD:0000070", "Stage 47"),
    ("CirobuD:0000071", "Stage 48"),
    ("CirobuD:0000072", "Stage 49"),
    ("CirobuD:0000073", "Stage 50"),
]

DICT_STAGES = {
    "CirobuD:0000021": 0,
    "CirobuD:0000022": 1,
    "CirobuD:0000023": 2,
    "CirobuD:0000024": 3,
    "CirobuD:0000025": 4,
    "CirobuD:0000026": 5,
    "CirobuD:0000027": 6,
    "CirobuD:0000028": 7,
    "CirobuD:0000029": 8,
    "CirobuD:0000030": 9,
    "CirobuD:0000031": 10,
    "CirobuD:0000032": 11,
    "CirobuD:0000033": 12,
    "CirobuD:0000034": 13,
    "CirobuD:0000035": 14,
    "CirobuD:0000036": 15,
    "CirobuD:0000037": 16,
    "CirobuD:0000038": 17,
    "CirobuD:0000039": 18,
    "CirobuD:0000040": 19,
    "CirobuD:0000041": 20,
    "CirobuD:0000042": 21,
    "CirobuD:0000043": 22,
    "CirobuD:0000044": 23,
    "CirobuD:0000045": 24,
    "CirobuD:0000046": 25,
    "CirobuD:0000047": 26,
    "CirobuD:0000048": 27,
    "CirobuD:0000049": 28,
    "CirobuD:0000050": 29,
    "CirobuD:0000051": 30,
    "CirobuD:0000052": 31,
    "CirobuD:0000053": 32,
    "CirobuD:0000054": 33,
    "CirobuD:0000055": 34,
    "CirobuD:0000056": 35,
    "CirobuD:0000057": 36,
    "CirobuD:0000058": 37,
    "CirobuD:0000059": 38,
    "CirobuD:0000060": 39,
    "CirobuD:0000061": 40,
    "CirobuD:0000062": 41,
    "CirobuD:0000063": 42,
    "CirobuD:0000064": 43,
    "CirobuD:0000070": 44,
    "CirobuD:0000071": 45,
    "CirobuD:0000072": 46,
    "CirobuD:0000073": 47,
}


class Term:
    def __init__(
        self,
        id: str,
        name: Optional[str] = None,
        namespace: Optional[str] = None,
        definition: Optional[str] = None,
        is_a: Optional[str] = None,
        relationship: Optional[Dict[str, List[str]]] = None,
        synonyms: Optional[List[str]] = None,
        references: Optional[List[str]] = None,
    ):
        self.id = id
        self.name = name
        self.namespace = namespace
        self.definition = definition
        self.is_a = is_a
        self.relationship = relationship or {}
        self.synonyms = synonyms or []
        self.references = references or []

    def __repr__(self):
        return (
            f"Term(id={self.id!r}, name={self.name!r}, namespace={self.namespace!r}, "
            f"definition={self.definition!r}, is_a={self.is_a!r}, "
            f"relationship={self.relationship!r}, synonyms={self.synonyms!r}, "
            f"references={self.references!r})"
        )


def parse_obo(
    filename: str,
) -> Tuple[List[Term], List[Tuple[str, str]], List[Tuple[str, str]]]:
    """
    Parses an OBO file to collect terms and their properties.

    Parameters:
    - filename (str): Path to the OBO file.

    Returns:
    - Tuple[List[Term], List[Tuple[str, str]], List[Tuple[str, str]]]: A tuple containing three items:
        1. A list of Term objects.
        2. A list of tuples representing part_of relationships where each tuple contains a term ID and the ID it is part of.
        3. A list of tuples representing preceded_by relationships where each tuple contains a term ID and the ID it is preceded by.
    """
    terms = []
    part_of = []
    preceded_by = []
    in_term = False
    attributes = {}

    with open(filename, "r") as f:
        for line in f:
            # Remove comment and strip whitespace
            value_comment = line.split("!", 1)
            line = value_comment[0].strip()

            if line == "[Term]":
                in_term = True
                attributes = {}
            elif line == "" and in_term:
                term = Term(
                    id=attributes.get("id"),
                    name=attributes.get("name"),
                    namespace=attributes.get("namespace"),
                    definition=attributes.get("def"),
                    is_a=attributes.get("is_a"),
                    relationship=attributes.get("relationship"),
                    synonyms=attributes.get("synonym"),
                    references=attributes.get("reference"),
                )
                terms.append(term)

                # Store part_of and preceded_by relationships separately as tuples
                if (
                    "relationship" in attributes
                    and "part_of" in attributes["relationship"]
                ):
                    for target_id in attributes["relationship"]["part_of"]:
                        part_of.append((term.id, target_id))
                if (
                    "relationship" in attributes
                    and "preceded_by" in attributes["relationship"]
                ):
                    for target_id in attributes["relationship"]["preceded_by"]:
                        preceded_by.append((term.id, target_id))

                in_term = False
                attributes = {}
            elif in_term:
                # Split at the first occurrence of ":"
                parts = line.split(":", 1)
                if len(parts) == 2:
                    tag = parts[0].strip()
                    value = parts[1].strip()

                    # Store relationships in a dictionary
                    if tag == "relationship":
                        rel_type, rel_value = value.split(" ", 1)
                        attributes.setdefault("relationship", {}).setdefault(
                            rel_type.strip(), []
                        ).append(rel_value.strip())
                    elif tag == "synonym":
                        attributes.setdefault("synonym", []).append(value)
                    elif tag == "reference":
                        attributes.setdefault("reference", []).append(value)
                    else:
                        attributes[tag] = value

    return terms, part_of, preceded_by


def split_hierarchy(terms: List[Term]) -> Tuple[List[Term], List[Term], List[Term]]:
    """
    Split a list of Term objects based on their names and separates them into three lists based on term ID prefixes.

    Parameters:
    - terms (List[Term]): The list of Term objects to be sorted.

    Returns:
    - Tuple[List[Term], List[Term], List[Term]]: A tuple containing three lists of Term objects, where the first list contains terms with ID starting with "CirobuA",
      the second list contains terms with ID starting with "CirobuD", and the third list contains the remaining terms.
    """
    cirobua_terms = []
    cirobud_terms = []
    other_terms = []

    for term in terms:
        if term.id.startswith("CirobuA"):
            cirobua_terms.append(term)
        elif term.id.startswith("CirobuD"):
            cirobud_terms.append(term)
        else:
            other_terms.append(term)

    return cirobua_terms, cirobud_terms, other_terms


def build_graph(filename: str) -> Tuple[nx.DiGraph, nx.DiGraph]:
    """
    Reads an OBO file and builds directed graphs based on the part_of and preceded_by relationships using networkx.
    Separates the terms into two different trees based on their ID prefixes.

    Parameters:
    - filename (str): Path to the OBO file.

    Returns:
    - Tuple[nx.DiGraph, nx.DiGraph]: Two directed graphs representing the OBO file based on the part_of and preceded_by relationships,
      where the first graph contains the CirobuA terms and the second graph contains the CirobuD terms.
    """
    terms, part_of, preceded_by = parse_obo(filename)

    G_cirobua = nx.DiGraph()
    G_cirobud = nx.DiGraph()

    cirobua_terms, cirobud_terms, _ = split_hierarchy(terms)

    for term in cirobua_terms:
        G_cirobua.add_node(term.id, name=term.name)
    for term in cirobud_terms:
        G_cirobud.add_node(term.id, name=term.name)

    for source, target in part_of:
        G_cirobua.add_edge(source, target, relationship="part_of")
    for source, target in preceded_by:
        G_cirobud.add_edge(source, target, relationship="preceded_by")

    return G_cirobua, G_cirobud


def get_linked_lists(graph: nx.DiGraph) -> List[List[Tuple[str, str]]]:
    """
    Returns a list of linked lists in the given directed graph.

    Parameters:
    - graph (nx.DiGraph): The directed graph.

    Returns:
    - List[List[Tuple[str, str]]]: A list of linked lists, where each linked list is represented as a list of tuples containing the node id and name in the order of appearance.
    """

    tree = nx.reverse(graph, copy=True)

    linked_lists = []
    visited = set()

    for node in tree.nodes:
        if node not in visited:
            linked_list = []
            current = node

            while current is not None:
                linked_list.append((current, tree.nodes[current]["name"]))
                visited.add(current)

                successors = list(tree.successors(current))
                if len(successors) == 1:
                    current = successors[0]
                else:
                    current = None

            linked_lists.append(linked_list)

    return linked_lists


def find_leaves(tree: nx.DiGraph, node: str) -> List[str]:
    """
    Recursively search for all leaves under a specific node in a tree.

    Parameters:
    - tree (nx.DiGraph): The tree to search.
    - node (str): The starting node.

    Returns:
    - List[str]: A list of all leaves found under the given node.
    """
    # If this node has no successors, it is a leaf.
    if tree.out_degree(node) == 0:
        return [node]

    leaves = []
    for successor in tree.successors(node):
        leaves.extend(find_leaves(tree, successor))

    return leaves


def node_to_leaves(graph: nx.DiGraph, node_attr="name") -> Dict[str, List[str]]:
    """
    Constructs a dictionary mapping from nodes to the list of leaves under them in the tree.
    Assumes that the graph represents a tree structure.

    Parameters:
    - graph (nx.DiGraph): The directed graph that needs to be processed.
    - node_attr (str): The node attribute used to fetch the node's name.

    Returns:
    - Dict[str, List[str]]: A dictionary mapping each node to a list of its leaves.
    """
    # First, we reverse the graph to make it a tree (parent -> children).
    tree = nx.reverse(graph, copy=True)
    # Compile regex pattern for efficiency and specificity
    pattern = re.compile("^([ABab])(\d+)\.(\d+)")

    node_to_leaves = {}
    for node in tree.nodes:
        if pattern.match(node_name := tree.nodes[node].get(node_attr)):
            leaves = find_leaves(tree, node)  # Find all leaves under this node.
            leaves_names = [
                name
                for leaf in leaves
                if pattern.match(name := tree.nodes[leaf].get(node_attr))
            ]  # Convert node IDs to names or use the node ID if name is not present.
            node_to_leaves[node_name] = leaves_names

    return node_to_leaves


def write_to_csv(data: Dict[str, List[str]], filename: str) -> None:
    """
    Writes the provided dictionary to a CSV file.

    Parameters:
    - data (Dict[str, List[str]]): The dictionary to write. Keys are strings, and values are lists of strings.
    - filename (str): The name of the output CSV file.
    """
    with open(filename, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)

        # Write header
        writer.writerow(["Node", "Leaves"])

        # Write data
        for node, leaves in data.items():
            writer.writerow([node, ",".join(leaves)])
