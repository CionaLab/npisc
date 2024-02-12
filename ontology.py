from typing import List, Tuple, Optional, Dict
import re
import csv
import networkx as nx


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


def build_graph(filename: str) -> nx.DiGraph:
    """
    Reads an OBO file and builds a directed graph based on the part_of relationship using networkx.

    Parameters:
    - filename (str): Path to the OBO file.

    Returns:
    - nx.DiGraph: A directed graph representation of the OBO file based on the part_of relationship.
    """
    terms, relationships = parse_obo(filename)

    G = nx.DiGraph()
    for term in terms:
        G.add_node(term.id, name=term.name)

    for source, target in relationships:
        G.add_edge(source, target, relationship="part_of")

    return G


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
