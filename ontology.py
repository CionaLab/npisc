from typing import List, Tuple, Optional, Dict, Callable
from collections import defaultdict
import csv
from urllib.parse import urlparse, parse_qs

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

    tuples_attr = defaultdict(list)

    in_term = False
    attributes = {}

    def parse_relationship(
        term: Term, attributes: Dict[str, Dict[str, str]], attr: str
    ) -> None:
        if "relationship" in attributes and attr in attributes["relationship"]:
            for target_id in attributes["relationship"][attr]:
                tuples_attr[attr].append((term.id, target_id))

    def clear_stack() -> None:

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
        parse_relationship(term, attributes, "part_of")
        parse_relationship(term, attributes, "preceded_by")

    with open(filename, "r") as f:
        for line in f:
            # Remove comment and strip whitespace
            value_comment = line.split("!", 1)
            line = value_comment[0].strip()

            if line == "[Term]":
                if in_term:
                    clear_stack()

                in_term = True
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

    if not attributes:
        clear_stack()

    return (
        terms,
        tuples_attr["part_of"],
        tuples_attr["preceded_by"],
    )


def split_terms(
    terms: List[Term], test_func: Callable[[Term], bool]
) -> Tuple[List[Term], List[Term]]:
    """
    Split a list of Term objects based on a test function and separates them into two lists.

    Parameters:
    - terms (List[Term]): The list of Term objects to be sorted.
    - test_func (Callable[[Term], bool]): The test function that takes a Term object as input and returns a boolean value.

    Returns:
    - Tuple[List[Term], List[Term]]: A tuple containing two lists of Term objects, where the first list contains terms for which the test function returns True,
      and the second list contains terms for which the test function returns False.
    """
    true_terms = []
    false_terms = []

    for term in terms:
        if test_func(term):
            true_terms.append(term)
        else:
            false_terms.append(term)

    return true_terms, false_terms


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

    cirobua_terms, rest = split_terms(terms, lambda x: x.id.startswith("CirobuA"))
    cirobud_terms, _ = split_terms(rest, lambda x: x.id.startswith("CirobuD"))

    for term in cirobua_terms:
        G_cirobua.add_node(term.id, name=term.name)
        try:
            G_cirobua.nodes[term.id]["start"] = term.relationship.get("start")[0]
            G_cirobua.nodes[term.id]["end"] = term.relationship.get("end")[0]
        except (KeyError, TypeError):
            pass
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


def node_to_leaves(graph: nx.DiGraph) -> Dict[str, List[str]]:
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

    node_to_leaves = {}
    for node in tree.nodes:
        leaves = find_leaves(tree, node)  # Find all leaves under this node.
        if leaves:
            node_to_leaves[node] = leaves
    return node_to_leaves


def group_leaves_by_time(
    graph: nx.DiGraph, leaves: Dict[str, List[str]]
) -> Dict[Tuple[str, str, str], List[str]]:
    """
    Groups the leaves in the given dictionary by their start and end nodes.

    Parameters:
    - graph (nx.DiGraph): The directed graph representing the ontology.
    - leaves (Dict[str, List[str]]): A dictionary mapping each node to a list of its leaves.

    Returns:
    - Dict[Tuple[str, str, str], List[str]]: A dictionary with (id, start, end) as key and a list of leaf ids as value.
    """

    results = defaultdict(list)
    starts = nx.get_node_attributes(graph, "start")
    ends = nx.get_node_attributes(graph, "end")

    for node, leaves in leaves.items():
        for leaf in leaves:
            start = starts.get(leaf)
            end = ends.get(leaf)
            results[(node, start, end)].append(leaf)

    return results


def get_node_name(graph: nx.DiGraph, node: str) -> Optional[str]:
    """
    Finds the name of a node in the graph based on its ID.

    Parameters:
    - graph (nx.DiGraph): The directed graph.
    - node (str): The ID of the node to find.

    Returns:
    - Optional[str]: The name of the node if found, None otherwise.
    """
    return graph.nodes[node].get("name")


def parse_url(url: str) -> dict:
    """
    Parses the given URL and returns a dictionary of its parameters.

    Parameters:
    - url (str): The URL to parse.

    Returns:
    - dict: A dictionary containing the parsed parameters.
    """
    parsed_url = urlparse(url)
    query_params = parse_qs(parsed_url.query)
    return {key: value[0] for key, value in query_params.items()}


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
