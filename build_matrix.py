"""
This module contains functions for analyzing gene expression data.
"""

from typing import Dict, Union, Tuple, Iterable, Callable
from itertools import tee
import re

from scipy.sparse import issparse
import pandas as pd
import anndata
from sklearn.metrics import pairwise_distances

from ontology import build_graph, node_to_leaves


def preprocess_tsv(
    filepath: str, mapping_filepath: str, obo_filepath: str
) -> pd.DataFrame:
    """
    Reads a TSV file and preprocesses it by extracting relevant information.
    It also maps the KH2012 gene model to the KY21 gene model using a provided mapping file.
    Uses ontology information from OBO file to map territories to underlying blastomeres.

    The function extracts the stage in parentheses, the KH number, and the expression territory
    (cell) from the TSV rows. Columns with NaN values are dropped.

    Parameters:
    - filepath (str): Path to the TSV file.
    - mapping_filepath (str): Path to the TSV file containing the mapping between KH2012 and KY21.
    - obo_filepath (str): Path to the OBO file.

    Returns:
    - pd.DataFrame: A preprocessed DataFrame.
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

    The DataFrame should have columns: Stage, Gene, and Territory.
    The returned adjacency matrix has cells in rows and genes in columns, indicating
    gene expression in the respective cells. If multiple edges exist between the same
    vertices in the adjacency list, they are considered as a single edge in the matrix.

    Parameters:
    - df (pd.DataFrame): A DataFrame containing the adjacency list.

    Returns:
    - pd.DataFrame: An adjacency matrix with cells in rows and genes in columns.
    """

    # Create the adjacency matrix using pivot_table
    # Fill NaN with 0 and convert float dtype to int for binary representation
    mat = df.pivot_table(
        index="Territory", columns="Gene", aggfunc="size", fill_value=0
    ).astype(int)

    # Convert any value greater than 0 to 1, to represent a single edge
    mat = (mat > 0).astype(int)

    return mat


def append_raw(adata: anndata.AnnData, adata_raw: anndata.AnnData) -> anndata.AnnData:
    """
    Append the observation data to raw data.

    Parameters:
    - adata (anndata.AnnData): The AnnData object to append the raw data to.
    - adata_raw (anndata.AnnData): The AnnData object containing the raw data to append.

    Returns:
    - anndata.AnnData: The AnnData object with the raw data appended.
    """

    adata_raw.obs = adata_raw.obs.merge(
        adata.obs,
        how="left",
        left_index=True,
        right_index=True,
        suffixes=("_raw", None),
    )

    return adata_raw


def to_df(obj: Union[anndata.AnnData, pd.DataFrame]) -> pd.DataFrame:
    """
    Convert an AnnData or DataFrame to a DataFrame.

    Parameters:
    - obj (Union[anndata.AnnData, pd.DataFrame]): The object to convert.

    Returns:
    - pd.DataFrame: The converted DataFrame.
    """
    if isinstance(obj, anndata.AnnData):
        if issparse(obj.X):
            return pd.DataFrame(
                obj.X.toarray(), columns=obj.var_names, index=obj.obs_names
            )
        return pd.DataFrame(obj.X, columns=obj.var_names, index=obj.obs_names)
    return obj


def to_obj(
    obj: Union[anndata.AnnData, pd.DataFrame], df: pd.DataFrame
) -> Union[anndata.AnnData, pd.DataFrame]:
    """
    Convert a DataFrame back to its original type (either AnnData or DataFrame).

    Parameters:
    - obj (Union[anndata.AnnData, pd.DataFrame]): The original object.
    - df (pd.DataFrame): The DataFrame to convert.

    Returns:
    - Union[anndata.AnnData, pd.DataFrame]: The converted object.
    """
    if isinstance(obj, anndata.AnnData):
        return anndata.AnnData(
            X=df.values, var=pd.DataFrame(index=df.columns), obs=obj.obs, uns=obj.uns
        )
    if isinstance(obj, pd.DataFrame):
        return df
    return None


def pad_compatible(
    obj1: Union[anndata.AnnData, pd.DataFrame],
    obj2: Union[anndata.AnnData, pd.DataFrame],
) -> Tuple[Union[anndata.AnnData, pd.DataFrame], Union[anndata.AnnData, pd.DataFrame]]:
    """
    Adjust both input objects (either AnnData or DataFrame) to have the same columns.

    Parameters:
    - obj1 (Union[anndata.AnnData, pd.DataFrame]): First object
    (either AnnData or DataFrame).
    - obj2 (Union[anndata.AnnData, pd.DataFrame]): Second object
    (either AnnData or DataFrame).

    Returns:
    - tuple: Tuple containing the adjusted objects.
    """

    df1 = to_df(obj1)
    df2 = to_df(obj2)

    # Get union of columns from both dataframes
    all_cols = df1.columns.union(df2.columns)

    # Reindex dataframes with all_cols and fill missing values with zeros
    df1 = df1.reindex(columns=all_cols, fill_value=0)
    df2 = df2.reindex(columns=all_cols, fill_value=0)

    # Now you can use this function in your code like this:
    obj1 = to_obj(obj1, df1)
    obj2 = to_obj(obj2, df2)

    return obj1, obj2


def split_adata(adata: anndata.AnnData, col: str) -> Dict[str, anndata.AnnData]:
    """
    Splits an AnnData object into a dictionary of AnnData objects based on a
    column in adata.obs.

    Parameters:
    - adata (anndata.AnnData): Input AnnData object.
    - col (str): The column name in adata.obs based on which the splitting
    should be done.

    Returns:
    - Dict[str, anndata.AnnData]: Dictionary with unique values from the column
    as keys and respective sub-AnnData objects as values.
    """

    return {value: adata[adata.obs[col] == value] for value in adata.obs[col].unique()}


def get_distance(
    obj1: Union[anndata.AnnData, pd.DataFrame],
    obj2: Union[anndata.AnnData, pd.DataFrame],
    metric: str | Callable = "euclidean",
) -> pd.DataFrame:
    """
    Compute the distance between two objects (either AnnData or DataFrame) using
    the specified metric.

    Parameters:
    - obj1 (Union[anndata.AnnData, pd.DataFrame]): The first object
    (either AnnData or DataFrame).
    - obj2 (Union[anndata.AnnData, pd.DataFrame]): The second object
    (either AnnData or DataFrame).
    - metric (str, optional): The distance metric to use. Defaults to "euclidean".

    Returns:
    - pd.DataFrame: A dataframe of distance values between the two objects.
    """

    df1 = to_df(obj1)
    df2 = to_df(obj2)

    # Ensure both dataframes have the same columns
    df1, df2 = pad_compatible(df1, df2)

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

    Parameters:
    - iterable (Iterable): Input iterable.
    - n (int): The number of elements in each tuple.

    Returns:
    - Iterable[Tuple[Iterable, ...]]: Iterable producing n-tuples of adjacent
    elements.
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

    Parameters:
    - territory (str): The territory string, e.g., "A9.32".

    Returns:
    - Tuple[str, int, str]: A tuple with cell lineage, rounds of division, and
    cell number.
    """
    match = re.match(r"^([ABab])(\d+)\.(\d+\**)$", territory)
    if match:
        lineage, rounds, cell_num = match.groups()
        return lineage, int(rounds), cell_num
    return None, None, None
