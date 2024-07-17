from typing import Dict, Union, Tuple, Iterable, Callable, List
from itertools import tee, product
import re
import numpy as np
from scipy.spatial import distance
from scipy.sparse import issparse
from scipy.stats import wasserstein_distance_nd
import pandas as pd
import anndata
from sklearn.preprocessing import normalize
from sklearn.metrics import pairwise_distances
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
    elif isinstance(obj, pd.DataFrame):
        return df


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
    Splits an AnnData object into a dictionary of AnnData objects based on a column in adata.obs.

    Parameters:
    - adata (anndata.AnnData): Input AnnData object.
    - col (str): The column name in adata.obs based on which the splitting should be done.

    Returns:
    - Dict[str, anndata.AnnData]: Dictionary with unique values from the column as keys and
                                  respective sub-AnnData objects as values.
    """

    return {value: adata[adata.obs[col] == value] for value in adata.obs[col].unique()}


def get_distance(
    obj1: Union[anndata.AnnData, pd.DataFrame],
    obj2: Union[anndata.AnnData, pd.DataFrame],
    metric: str | Callable = "euclidean",
) -> pd.DataFrame:
    """
    Compute the distance between two objects (either AnnData or DataFrame) using the specified metric.

    Parameters:
    - obj1 (Union[anndata.AnnData, pd.DataFrame]): The first object (either AnnData or DataFrame).
    - obj2 (Union[anndata.AnnData, pd.DataFrame]): The second object (either AnnData or DataFrame).
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


def get_mahalanobis_distance(
    obj1: Union[anndata.AnnData, pd.DataFrame],
    obj2: Union[anndata.AnnData, pd.DataFrame],
    output_similarity: bool = False,
) -> pd.DataFrame:
    """
    Compute the Mahalanobis distance between two objects (either AnnData or DataFrame).
    Each RNAseq cell is a row vector. Covariances are calculated from obj1.

    Parameters:
    - obj1 (Union[anndata.AnnData, pd.DataFrame]): The first object (either AnnData or DataFrame).
    - obj2 (Union[anndata.AnnData, pd.DataFrame]): The second object (either AnnData or DataFrame).
    - output_similarity (bool, optional): If True, output the similarity matrix instead of the distance matrix.
                                          Defaults to False.

    Returns:
    - pd.DataFrame: A dataframe of Mahalanobis distance or similarity values.
    """

    df1 = to_df(obj1)
    df2 = to_df(obj2)

    # Ensure both dataframes have the same columns
    df1, df2 = pad_compatible(df1, df2)

    # Convert dataframes to numpy arrays
    arr1 = df1.to_numpy()
    arr2 = df2.to_numpy()

    # Calculate the gene experience covariance matrix
    cov = np.cov(arr1, rowvar=False)
    cov_i = np.linalg.pinv(cov)

    df_distance = pd.DataFrame(
        pairwise_distances(arr1, arr2, metric="mahalanobis", VI=cov_i, n_jobs=-1),
        index=df1.index,
        columns=df2.index,
    )

    df_distance = df_distance.sort_index(axis=0).sort_index(axis=1)

    if not output_similarity:
        return df_distance
    else:
        return 1 - (df_distance / df_distance.to_numpy().max())


def adjacent(iterable: Iterable, n: int = 2) -> Iterable[Tuple[Iterable, ...]]:
    """
    Return n-tuples of adjacent elements from the input iterable.

    Parameters:
    - iterable (Iterable): Input iterable.
    - n (int): The number of elements in each tuple.

    Returns:
    - Iterable[Tuple[Iterable, ...]]: Iterable producing n-tuples of adjacent elements.
    """
    iterators = tee(iterable, n)
    for i, iterator in enumerate(iterators):
        for _ in range(i):
            next(iterator, None)
    return zip(*iterators)


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
        (round1, round2): 1
        - get_distance(adj_matrices[round1], adj_matrices[round2], metric="cosine")
        for round1, round2 in adjacent(grouped.keys())
    }

    return similarities


def find_similar_clusters(
    adata1: anndata.AnnData, adata2: anndata.AnnData
) -> pd.DataFrame:
    """
    Find the most similar Leiden cluster pairs between two AnnData objects.

    Parameters:
    - adata1 (anndata.AnnData): The first AnnData object with Leiden clustering.
    - adata2 (anndata.AnnData): The second AnnData object with Leiden clustering.

    Returns:
    - pd.DataFrame: A dataframe with the most similar Leiden cluster pairs and their Wasserstein distance.
    """

    adata1, adata2 = pad_compatible(adata1, adata2)
    split_adata1 = split_adata(adata1, col="leiden")
    split_adata2 = split_adata(adata2, col="leiden")

    split_adata1 = {name: to_df(data) for name, data in split_adata1.items()}
    split_adata2 = {name: to_df(data) for name, data in split_adata2.items()}

    distances = [
        (name1, name2, wasserstein_distance_nd(data1, data2))
        for (name1, data1), (name2, data2) in product(
            split_adata1.items(), split_adata2.items()
        )
    ]

    df = pd.DataFrame(distances, columns=["leiden_1", "leiden_2", "distance"])

    return df.sort_values(by="distance").reset_index(drop=True)


def find_coi(
    df: pd.DataFrame,
    adata: anndata.AnnData,
    cluster: str = "leiden",
    blastomere: str = "Territory_eq",
    quantile: float = 0.4,
) -> List:
    """
    Find the cells of interest based on the expression similarity to blastomeres in the DataFrame.

    Parameters:
    - df (pandas.DataFrame): The input DataFrame containing the similarity matrix.
    - adata (anndata.AnnData): The input AnnData object.
    - cluster (str, optional): The column name in adata.obs based on which the clusters of interest are determined. Default is "leiden".
    - blastomere (str, optional): The column name in df representing the blastomere. Default is "Territory_eq".
    - quantile (float, optional): The quantile value used to determine the similarity threshold. Default is 0.4.

    Returns:
    - List: A list of cluster labels that have expression similarity to blastomeres above the specified quantile threshold.

    """
    df = (
        df.melt(ignore_index=False)
        .join(adata.obs[cluster])
        .groupby(blastomere)
        .apply(
            lambda x: x[x["value"] <= x["value"].quantile(quantile)],
            include_groups=False,
        )
        .reset_index(level=blastomere)
    )
    return df.index.astype(str).unique().to_list()


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
