"""
This module contains functions for analyzing gene expression data.
"""

from typing import Dict, List, Tuple, Iterable, Callable
from itertools import tee
import re

import numpy as np
import pandas as pd
import anndata
from sklearn.metrics import pairwise_distances
import matplotlib.axes
import geopandas as gpd
import networkx as nx

from .ontology import build_graph, node_to_leaves


def preprocess_tsv(
    filepath: str, mapping_filepath: str, obo_filepath: str
) -> pd.DataFrame:
    """
    Reads a TSV file and preprocesses it by extracting relevant information.
    It also maps the KH2012 gene model to the KY21 gene model using a provided
    mapping file.

    Uses ontology information from OBO file to map territories to underlying
    blastomeres.

    The function extracts the stage in parentheses, the KH number, and the
    expression territory (cell) from the TSV rows. Columns with NaN values are
    dropped.

    :param filepath: Path to the TSV file.
    :type filepath: str
    :param mapping_filepath: Path to the TSV file containing the mapping between
        KH2012 and KY21.
    :type mapping_filepath: str
    :param obo_filepath: Path to the OBO file.
    :type obo_filepath: str
    :return: A preprocessed DataFrame.
    :rtype: pd.DataFrame
    """

    # Build graph from OBO file
    graph = build_graph(obo_filepath)
    # Create the mapping from nodes to their leaf nodes
    map_blastomeres = node_to_leaves(graph)

    # Read the TSV file into a DataFrame
    df = pd.read_csv(
        filepath, sep="\t", header=None, names=["URL", "Stage", "Gene", "Territory"]
    )

    # Extract stage in parentheses
    df["Stage"] = df["Stage"].str.extract(r"\((.*?)\)", expand=False)

    # Extract the KH number using regex
    df["Gene"] = df["Gene"].str.extract(r"(KH2012:KH\.[A-Z]\d+\.\d+)", expand=False)

    # Convert the Territory column to corresponding blastomeres using the ontology information
    df["Territory"] = df["Territory"].map(map_blastomeres)
    df = df.explode(column="Territory")

    # Load the mapping TSV into a DataFrame
    mapping_df = pd.read_csv(
        mapping_filepath, sep="\t", header=None, names=["KH2012", "KY21"]
    )

    # Add the prefixes back to the KH2012 and KY21 columns
    mapping_df["KH2012"] = "KH2012:" + mapping_df["KH2012"]
    mapping_df["KY21"] = "KY21:" + mapping_df["KY21"]

    # Merge the original DataFrame with the mapping DataFrame to map the genes
    df = pd.merge(df, mapping_df, how="left", left_on="Gene", right_on="KH2012")

    # Replace the Gene column with the KY21 mapping
    df["Gene"] = df["KY21"]

    # Drop rows with NaN values and the extra columns from the merge
    df = df[["URL", "Stage", "Gene", "Territory"]].dropna()

    return df


def build_from_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    Converts a DataFrame containing an adjacency list into an adjacency matrix.

    The DataFrame should have columns: Stage, Gene, and Territory.  The returned
    adjacency matrix has cells in rows and genes in columns, indicating gene
    expression in the respective cells. If multiple edges exist between the same
    vertices in the adjacency list, they are considered as a single edge in the
    matrix.

    :param df: A DataFrame containing the adjacency list.
    :type df: pd.DataFrame
    :return: An adjacency matrix with cells in rows and genes in columns.
    :rtype: pd.DataFrame
    """

    # Create the adjacency matrix using pivot_table
    # Fill NaN with 0 and convert float dtype to int for binary representation
    mat = df.pivot_table(
        index="Territory", columns="Gene", aggfunc="size", fill_value=0
    ).astype(int)

    # Convert any value greater than 0 to 1, to represent a single edge
    mat = (mat > 0).astype(int)

    return mat


def pad_df_adata(
    df: pd.DataFrame,
    adata: anndata.AnnData,
) -> pd.DataFrame:
    """
    Adjust the DataFrame to have the same columns by padding the missing columns
    with zeros and removing the extra columns.

    :param df: The dataframe to be adjusted.
    :type df: pd.DataFrame
    :param adata: The reference AnnData object.
    :type df2: annadata.AnnData
    :return: the adjusted DataFrames.
    :rtype: pd.DataFrame
    """

    # Reindex dataframes with all_cols and fill missing values with zeros
    df = df.reindex(columns=adata.var_names, fill_value=0)

    return df


def pad_dfs(
    df1: pd.DataFrame,
    df2: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Adjust both DataFrame to have the same
    columns.

    :param df1: First DataFrame.
    :type df1: pd.DataFrame
    :param df2: Second DataFrame.
    :type df2: pd.DataFrame
    :return: Tuple containing the adjusted DataFrames.
    :rtype: Tuple[pd.DataFrame, pd.DataFrame]
    """

    # Get union of columns from both dataframes
    all_cols = df1.columns.union(df2.columns)

    # Reindex dataframes with all_cols and fill missing values with zeros
    df1 = df1.reindex(columns=all_cols, fill_value=0)
    df2 = df2.reindex(columns=all_cols, fill_value=0)

    return df1, df2


def split_adata(adata: anndata.AnnData, col: str) -> Dict[str, anndata.AnnData]:
    """
    Splits an AnnData object into a dictionary of AnnData objects based on a
    column in adata.obs.

    :param adata: Input AnnData object.
    :type adata: anndata.AnnData
    :param col: The column name in adata.obs based on which the splitting should
        be done.
    :type col: str
    :return: Dictionary with unique values from the column as keys and
        respective sub-AnnData objects as values.
    :rtype: Dict[str, anndata.AnnData]
    """

    return {value: adata[adata.obs[col] == value] for value in adata.obs[col].unique()}


def get_distance(
    df1: pd.DataFrame,
    df2: pd.DataFrame,
    metric: str | Callable = "euclidean",
) -> pd.DataFrame:
    """
    Compute the distance between two DataFrames using the specified metric.

    :param df1: The first DataFrames.
    :type df1: pd.DataFrame
    :param df2: The second DataFrame.
    :type df2: pd.DataFrame
    :param metric: The distance metric to use. Defaults to "euclidean".
    :type metric: str, optional
    :return: A dataframe of distance values between the two objects.
    :rtype: pd.DataFrame
    """

    # Convert dataframes to numpy arrays
    arr1 = df1.to_numpy()
    arr2 = df2.to_numpy()

    similarity_df = pd.DataFrame(
        pairwise_distances(arr1, arr2, metric=metric),
        index=df1.index,
        columns=df2.index,
    )

    # Sort the result
    similarity_df = similarity_df.sort_index(axis=0).sort_index(axis=1)

    return similarity_df


def adjacent(iterable: Iterable, n: int = 2) -> Iterable[Tuple[Iterable, ...]]:
    """
    Return n-tuples of adjacent elements from the input iterable.

    :param iterable: Input iterable.
    :type iterable: Iterable
    :param n: The number of elements in each tuple.
    :type n: int
    :return: Iterable producing n-tuples of adjacent elements.
    :rtype: Iterable[Tuple[Iterable, ...]]
    """

    iterators = tee(iterable, n)
    for i, iterator in enumerate(iterators):
        for _ in range(i):
            next(iterator, None)
    return zip(*iterators)


def parse_territory(territory: str) -> Tuple[str, int, str]:
    """
    Parse the territory string into its components: cell lineage, rounds of
    division, and cell number.

    :param territory: The territory string, e.g., "A9.32".
    :type territory: str
    :return: A tuple with cell lineage, rounds of division, and cell number.
    :rtype: tuple[str, int, str]
    """

    match = re.match(r"^([ABab])(\d+)\.(\d+\**)$", territory)
    if match:
        lineage, rounds, cell_num = match.groups()
        return lineage, int(rounds), cell_num
    return None, None, None


def map_cells(
    adata: anndata.AnnData,
    pattern: np.ndarray,
    basis: str = "PCs",
    use_rep: str = "X_pca",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Maps scRNAseq cells to blastomeres using a given pattern.

    :param adata: Annotated data matrix.
    :type adata: anndata.AnnData
    :param pattern: Pattern to be used for mapping.
    :type pattern: np.ndarray
    :param basis: Basis to be used from adata.varm. Default is "PCs".
    :type basis: str
    :param use_rep: Representation to be used from adata.obsm. Default is "X_pca".
    :type use_rep: str
    :return: A tuple with a dataframe with the cosine similarity of individual
        scRNAseq cells, a dataframe with the mean cosine similarity of each
        Leiden cluster, and the most similar blastomere for each Leiden cluster.
    :rtype: tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]
    """

    pattern = pad_df_adata(pattern, adata)

    df_cells = 1 - get_distance(
        pattern @ adata.varm[basis],
        pd.DataFrame(adata.obsm[use_rep], index=adata.obs_names),
        "cosine",
    )

    df_leiden = (
        df_cells.T.join(adata.obs["leiden"], how="left").groupby("leiden").mean()
    )

    df_map = pd.DataFrame(
        {"cluster": df_leiden.idxmax(axis=0), "cos_theta": df_leiden.max(axis=0)}
    )

    return df_cells, df_leiden, df_map


def plot_np(
    df: pd.DataFrame,
    gdf: gpd.GeoDataFrame,
    ax: matplotlib.axes.Axes,
    gpd_kwds: dict,
    anno_kwds: dict,
) -> matplotlib.axes.Axes:
    """
    Plot the most similar Leiden cluster on a blastomere map. The merged
    dataframe should have a cluster column for annotation.

    :param df: The DataFrame containing the Leiden cluster and the most similar
    blastomere.
    :type df: pd.DataFrame
    :param gdf: The GeoDataFrame containing the blastomere map.
    :type gdf: gpd.GeoDataFrame
    :param ax: The Axes for the plot.
    :type ax: matplotlib.axes.Axes
    :param gpd_kwds: Keyword arguments to pass to the GeoDataFrame plot function.
        It should contain the column name to use for coloring.
    :type gpd_kwds: dict
    :param anno_kwds: Keyword arguments to pass to the annotate function.
    :type anno_kwds: dict
    :return: The axes object representing the plot.
    :rtype: matplotlib.axes.Axes
    """

    gdf = gdf.merge(df, left_on="name", right_index=True, how="left")
    gdf.plot(ax=ax, **gpd_kwds)
    gdf.apply(
        lambda x: ax.annotate(
            text=f"{x['cluster']}",
            xy=x.geometry.centroid.coords[0],
            ha="center",
            **anno_kwds,
        ),
        axis=1,
    )
    ax.axis("off")

    return ax


def cross_stage_distance(
    adata1: anndata.AnnData,
    adata2: anndata.AnnData,
    name1: str,
    name2: str,
    use_rep: str = "X_pca",
) -> pd.DataFrame:
    """
    Calculate the cross-stage distance between two AnnData objects.

    :param adata1: The first AnnData object.
    :type adata1: AnnData
    :param adata2: The second AnnData object.
    :type adata2: AnnData
    :param name1: The name of the second AnnData object.
    :type name1: str
    :param name2: The name of the second AnnData object.
    :type name2: str
    :param use_rep: Representation to be used from adata.obsm. Default is "X_pca".
    :type use_rep: str
    :return: The round trip distance matrix.
    :rtype: pd.DataFrame
    """

    t1 = get_distance(
        pd.DataFrame(adata1.obsm[use_rep], index=adata1.obs_names),
        pd.DataFrame(adata2.obsm[use_rep], index=adata2.obs_names),
        "cosine",
    )

    d1 = (
        t1.melt(ignore_index=False, var_name="target", value_name="cos_theta")
        .reset_index(names="source")
        .merge(adata1.obs["leiden"], left_on=["source"], right_index=True)
        .merge(
            adata2.obs["leiden"],
            left_on=["target"],
            right_index=True,
            suffixes=(f"_{name1}", f"_{name2}"),
        )
        .groupby([f"leiden_{name1}", f"leiden_{name2}"])["cos_theta"]
        .mean()
        .reset_index()
    )

    d3 = d1.pivot(
        index=f"leiden_{name1}", columns=f"leiden_{name2}", values="cos_theta"
    )

    return d3


def make_stage_name(stage: str) -> Callable[[str], str]:
    """
    Creates a function that appends a given cell name to a specified stage.

    :param stage: The stage to be prefixed to the cell name.
    :type stage: str
    :return: A function that takes a cell name as input and returns a string
        combining the stage and the cell name in the format "{stage}_{cell}".
    :rtype: Callable[[str], str]
    """

    def make_name(cell: str) -> str:
        return f"{stage}_{cell}"

    return make_name


def make_digraph(
    d_adatas: Dict[str, anndata.AnnData],
    stages: Iterable,
    use_rep: str = None,
) -> nx.DiGraph:
    """
    Create a directed graph (DiGraph) from single-cell data across different
    stages.

    This function constructs a directed graph where nodes represent cells from
    different stages, and edges represent the distance between cells from
    adjacent stages.

    :param d_adatas: Dictionary of AnnData objects, where keys are stage
        identifiers and values are AnnData objects containing single-cell data
        for each stage.
    :type d_adatas: Dict[anndata.AnnData]
    :param stages: An iterable of stage identifiers, defining the order of
        stages.
    :type stages: Iterable
    :param use_rep: Optional. The representation to use for computing distances.
        If None, the default representation is used.
    :type use_rep: str, optional
    :return: A directed graph with weighted edges representing distances between
        cells from adjacent stages.
    :rtype: nx.DiGraph
    """

    dg = nx.DiGraph()

    for k1, k2 in adjacent(stages):

        df = cross_stage_distance(
            d_adatas[k1],
            d_adatas[k2],
            k1,
            k2,
            use_rep=use_rep,
        )

        edges = df.stack().reset_index()
        edges.columns = ["source", "target", "weight"]
        edges["source"] = edges["source"].apply(make_stage_name(k1))
        edges["target"] = edges["target"].apply(make_stage_name(k2))
        dg.add_weighted_edges_from(edges.values)

    return dg


class CrossStage:
    """
    CrossStage class for managing and updating shortest paths in a directed
    graph.

    :param graph: The directed graph on which shortest paths are calculated.
    :type graph: nx.DiGraph
    :param nodes_start: A list of starting node identifiers.
    :type nodes_start: List[str]
    :param nodes_end: A list of ending node identifiers.
    :type nodes_end: List[str]
    """

    def __init__(self, graph: nx.DiGraph, nodes_start: List[str], nodes_end: List[str]):
        self.graph = graph
        self.update_internal(nodes_start, nodes_end)

    def update_internal(self, nodes_start: List[str], nodes_end: List[str]):
        """
        Updates the internal shortest paths data structures for the given start
            and end nodes.

        This method calculates the shortest paths between each pair of nodes
        from `nodes_start` to `nodes_end`
        using Dijkstra's algorithm with edge weights. It updates two internal
        attributes:
        - ``self.df_shortest_paths``: A DataFrame containing the shortest path
            distances from each source node
        - ``self.d_shortest_paths``: A dictionary where keys are tuples of
            (source, target) nodes and values

        :param nodes_start: A list of starting node identifiers.
        :type nodes_start: List[str]
        :param nodes_end: A list of ending node identifiers.
        :type nodes_end: List[str]
        """

        self.df_shortest_paths = (
            pd.DataFrame(
                [
                    {
                        "source": s,
                        "target": t,
                        "distance": nx.shortest_path_length(
                            self.graph,
                            source=s,
                            target=t,
                            weight="weight",
                        ),
                    }
                    for s in nodes_start
                    for t in nodes_end
                ]
            )
            .groupby("source")
            .apply(lambda x: x.nsmallest(1, "distance"))
            .reset_index(drop=True)
            .set_index("source")
        )

        self.d_shortest_paths = {
            (s, t): nx.shortest_path(
                self.graph,
                source=s,
                target=t,
                weight="weight",
            )
            for s in nodes_start
            for t in nodes_end
        }

    def find_paths(self, source: str) -> List[str]:
        """
        Finds the shortest paths from a given source node to all target nodes.

        :param source: The source node identifier.
        :type source: str
        :return: A list containing the nodes on the shortest paths from the
            source node.
        :rtype: List[str]
        """
        return self.d_shortest_paths[
            source, self.df_shortest_paths.loc[source]["target"]
        ]
