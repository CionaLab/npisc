from typing import Dict, Union, Tuple, Iterable
from itertools import tee
import re
import numpy as np
from scipy.spatial import distance
import pandas as pd
import anndata
from sklearn.preprocessing import normalize
import seaborn as sns
import matplotlib.pyplot as plt

from .ontology import build_graph, node_to_leaves


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


def pad_compatible(
    obj1: Union[anndata.AnnData, pd.DataFrame],
    obj2: Union[anndata.AnnData, pd.DataFrame],
) -> Tuple[Union[anndata.AnnData, pd.DataFrame], Union[anndata.AnnData, pd.DataFrame]]:
    """
    Adjust both input objects (either AnnData or DataFrame) to have the same columns.

    Parameters:
    - obj1 (Union[anndata.AnnData, pd.DataFrame]): First object (either AnnData or DataFrame).
    - obj2 (Union[anndata.AnnData, pd.DataFrame]): Second object (either AnnData or DataFrame).

    Returns:
    - tuple: Tuple containing the adjusted objects.
    """

    # Convert AnnData to DataFrame if necessary
    df1 = (
        pd.DataFrame(obj1.X, columns=obj1.var_names, index=obj1.obs_names)
        if isinstance(obj1, anndata.AnnData)
        else obj1
    )
    df2 = (
        pd.DataFrame(obj2.X, columns=obj2.var_names, index=obj2.obs_names)
        if isinstance(obj2, anndata.AnnData)
        else obj2
    )

    # Get union of columns from both dataframes
    all_cols = df1.columns.union(df2.columns)

    # Reindex dataframes with all_cols and fill missing values with zeros
    df1 = df1.reindex(columns=all_cols, fill_value=0)
    df2 = df2.reindex(columns=all_cols, fill_value=0)

    # Convert back to AnnData if original objects were AnnData
    if isinstance(obj1, anndata.AnnData):
        obj1.X = df1.values
        obj1.var_names = all_cols
    elif isinstance(obj1, pd.DataFrame):
        obj1 = df1

    if isinstance(obj2, anndata.AnnData):
        obj2.X = df2.values
        obj2.var_names = all_cols
    elif isinstance(obj2, pd.DataFrame):
        obj2 = df2

    return obj1, obj2


def split_adata(adata: anndata.AnnData, col: str) -> Dict[str, anndata.AnnData]:
    """
    Splits an AnnData object into a dictionary of AnnData objects based on a column in adata.obs.

    Parameters:
    - adata (anndata.AnnData): Input AnnData object.
    - col (str): The column name in adata.obs based on which the splitting should be done.

    Returns:
    - Dict[str, anndata.AnnData]: Dictionary with unique values from the column as keys and
                                  respective sub-AnnData objects as values.
    """

    return {value: adata[adata.obs[col] == value] for value in adata.obs[col].unique()}


def get_cos_similarity(
    obj1: Union[anndata.AnnData, pd.DataFrame],
    obj2: Union[anndata.AnnData, pd.DataFrame],
) -> pd.DataFrame:
    """
    Compute the cosine similarity between two objects (either AnnData or DataFrame).

    Parameters:
    - obj1 (Union[anndata.AnnData, pd.DataFrame]): First object (either AnnData or DataFrame).
    - obj2 (Union[anndata.AnnData, pd.DataFrame]): Second object (either AnnData or DataFrame).

    Returns:
    - pd.DataFrame: A dataframe of cosine similarity values.
    """

    # Convert AnnData to DataFrame if necessary
    df1 = (
        pd.DataFrame(obj1.X, columns=obj1.var_names, index=obj1.obs_names)
        if isinstance(obj1, anndata.AnnData)
        else obj1
    )
    df2 = (
        pd.DataFrame(obj2.X, columns=obj2.var_names, index=obj2.obs_names)
        if isinstance(obj2, anndata.AnnData)
        else obj2
    )

    # Ensure both dataframes have the same columns
    df1, df2 = pad_compatible(df1, df2)

    # Convert dataframes to numpy arrays
    arr1 = df1.to_numpy()
    arr2 = df2.to_numpy()

    # Normalize both arrays using L2 norm
    arr1 = normalize(arr1, axis=1, norm="l2")
    arr2 = normalize(arr2, axis=1, norm="l2")

    # Compute cosine similarity
    similarity = arr1 @ arr2.T

    # Convert the result back to a dataframe with appropriate row and column names
    similarity_df = pd.DataFrame(similarity, index=df1.index, columns=df2.index)

    # Sort the result
    similarity_df = similarity_df.sort_index(axis=0).sort_index(axis=1)

    return similarity_df

def get_mahalanobis_distance(
    obj1: Union[anndata.AnnData, pd.DataFrame],
    obj2: Union[anndata.AnnData, pd.DataFrame],
) -> pd.DataFrame:
    """
    Compute the Mahalanobis distance between two objects (either AnnData or DataFrame).

    Parameters:
    - obj1 (Union[anndata.AnnData, pd.DataFrame]): First object (either AnnData or DataFrame).
    - obj2 (Union[anndata.AnnData, pd.DataFrame]): Second object (either AnnData or DataFrame).

    Returns:
    - pd.DataFrame: A dataframe of Mahalanobis distance values.
    """

    # Convert AnnData to DataFrame if necessary
    df1 = (
        pd.DataFrame(obj1.X, columns=obj1.var_names, index=obj1.obs_names)
        if isinstance(obj1, anndata.AnnData)
        else obj1
    )
    df2 = (
        pd.DataFrame(obj2.X, columns=obj2.var_names, index=obj2.obs_names)
        if isinstance(obj2, anndata.AnnData)
        else obj2
    )

    # Ensure both dataframes have the same columns
    df1, df2 = pad_compatible(df1, df2)

    # Convert dataframes to numpy arrays
    arr1 = df1.to_numpy()
    arr2 = df2.to_numpy()

    # Calculate the covariance matrix
    cov_matrix = np.cov(arr1.T)
    inv_cov_matrix = np.linalg.inv(cov_matrix)

    # Compute Mahalanobis distance
    distance_matrix = distance.cdist(arr1, arr2, metric='mahalanobis', VI=inv_cov_matrix)

    # Convert the result back to a dataframe with appropriate row and column names
    distance_df = pd.DataFrame(distance_matrix, index=df1.index, columns=df2.index)

    # Sort the result
    distance_df = distance_df.sort_index(axis=0).sort_index(axis=1)

    return distance_df


def pairwise(iterable: Iterable) -> Iterable[Tuple[int, int]]:
    """
    Return pairs of adjacent elements from the input iterable.

    Parameters:
    - iterable (Iterable): Input iterable.

    Returns:
    - Iterable[Tuple[int, int]]: Iterable producing pairs of adjacent elements.
    """
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def parse_territory(territory: str) -> Tuple[str, int, str]:
    """
    Parse the territory string into its components: cell lineage, rounds of division, and cell number.

    Parameters:
    - territory (str): The territory string, e.g., "A9.32".

    Returns:
    - Tuple[str, int, str]: A tuple with cell lineage, rounds of division, and cell number.
    """
    match = re.match(r"^([ABab])(\d+)\.(\d+\**)$", territory)
    if match:
        lineage, rounds, cell_num = match.groups()
        return lineage, int(rounds), cell_num
    else:
        return None, None, None


def compute_stage_similarity(df: pd.DataFrame) -> Dict[Tuple[int, int], pd.DataFrame]:
    """
    Compute the cosine similarity between adjacency matrices of two adjacent division rounds.

    Parameters:
    - df (pd.DataFrame): The dataframe returned by preprocess_tsv.

    Returns:
    - Dict[Tuple[int, int], pd.DataFrame]: A dictionary where keys are tuples with two division rounds,
                                           and values are dataframes of cosine similarity values between
                                           the adjacency matrices of the two rounds.
    """

    # Apply the parse_territory function to the "Territory" column
    df["Lineage"], df["Rounds"], df["Cell_Num"] = zip(
        *df["Territory"].apply(parse_territory)
    )

    # Group the dataframe by the 'Rounds' column
    grouped = {name: group for name, group in df.groupby("Rounds")}

    # Build the adjacency matrix for each group
    adj_matrices = {
        round_num: build_from_df(round_df) for round_num, round_df in grouped.items()
    }

    # Calculate cosine similarity for adjacent groups using pairwise and dictionary comprehension
    similarities = {
        (round1, round2): get_cos_similarity(adj_matrices[round1], adj_matrices[round2])
        for round1, round2 in pairwise(grouped.keys())
    }

    return similarities


def find_similar_clusters(
    adata1: anndata.AnnData, adata2: anndata.AnnData
) -> pd.DataFrame:
    """
    Find the most similar Leiden cluster pairs between two AnnData objects.

    Parameters:
    - adata1 (anndata.AnnData): First AnnData object with Leiden clustering.
    - adata2 (anndata.AnnData): Second AnnData object with Leiden clustering.

    Returns:
    - pd.DataFrame: A dataframe with the most similar Leiden cluster pairs and their average cosine similarity.
    """

    # Calculate cosine similarity using the existing function
    cos_sim = get_cos_similarity(adata1, adata2)

    # Group by Leiden clusters in adata1
    group1 = adata1.obs["leiden"]
    mean_sim1 = cos_sim.groupby(group1, axis=0).mean()

    # Group by Leiden clusters in adata2
    group2 = adata2.obs["leiden"]
    mean_sim2 = mean_sim1.groupby(group2, axis=1).mean()

    # Flatten the dataframe and reset index
    flattened = mean_sim2.reset_index().melt(
        id_vars="leiden", var_name="leiden_2", value_name="similarity"
    )

    # Sort by similarity and return
    return flattened.sort_values(by="similarity", ascending=False).reset_index(
        drop=True
    )


def plot_distance(matrix: pd.DataFrame) -> None:
    """
    Plot a heatmap and dendrogram of a distance matrix.

    Parameters:
    - matrix (pd.DataFrame): DataFrame with the distance matrix.

    Returns:
    - None: Displays a heatmap and dendrogram.
    """

    # Fill diagonal and NaN values (if any) with 0 for better visualization
    np.fill_diagonal(matrix.values, 0)
    matrix = matrix.fillna(0)

    # Plot heatmap and dendrogram using seaborn's clustermap
    sns.clustermap(matrix, cmap="viridis", xticklabels=True, yticklabels=True)

    # Display the plot
    plt.show()
