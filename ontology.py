from typing import List, Tuple, Optional, Dict
import csv
import networkx as nx


def parse_obo(filename: str) -> Tuple[Dict[str, str], List[Tuple[str, str]]]:
    """
    Parses an OBO file for nodes and relationships.

    Parameters:
    - filename (str): Path to the OBO file.

    Returns:
    - Tuple[Dict[str, str], List[Tuple[str, str]]]:
      - Dictionary of term IDs to term names.
      - List of (source, target) tuples representing part_of relationships.
    """
    terms = {}
    relationships = []

    in_term = False
    attributes = {}  # Dictionary to store parsed attributes for the current term

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if line == "[Term]":
                in_term = True
                attributes = {}  # Reset attributes dictionary for the new term
            elif line == "" and in_term:
                # Only consider nodes that match our regex pattern
                if "name" in attributes:
                    terms[attributes.get("id")] = attributes.get("name")

                # If there's a part_of relationship, add it to our relationships
                if (
                    "relationship" in attributes
                    and "part_of" in attributes["relationship"]
                ):
                    relationships.append(
                        (attributes["id"], attributes["relationship"].split()[1])
                    )

                in_term = False
            elif in_term:
                # Parse line based on the "tag: value ! comment" specification
                parts = line.split(": ", 1)  # Split at first occurrence of ": "
                if len(parts) == 2:
                    tag = parts[0].strip()
                    value_comment = parts[1].split(" ! ")
                    value = value_comment[0].strip()

                    # Depending on the tag, process and store the value
                    if tag == "id":
                        attributes["id"] = value
                    elif tag == "name":
                        attributes["name"] = value
                    elif tag == "relationship" and "part_of" in value:
                        attributes["relationship"] = value

    return terms, relationships


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
    for term_id, term_name in terms.items():
        G.add_node(term_id, name=term_name)

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

    node_to_leaves = {}
    for node in tree.nodes:
        node_name = tree.nodes[node].get(
            node_attr, node
        )  # Get the name of the node or use node ID if name isn't present.
        leaves = find_leaves(tree, node)  # Find all leaves under this node.
        leaves_names = [
            tree.nodes[leaf].get(node_attr, leaf) for leaf in leaves
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
