from typing import Dict
import numpy as np
import pandas as pd
import anndata


def preprocess_tsv(filepath: str, mapping_filepath: str) -> pd.DataFrame:
    """
    Reads a TSV file and preprocesses it by extracting relevant information.
    It also maps the KH2012 gene model to the KY21 gene model using a provided mapping file.

    The function extracts the stage in parentheses, the KH number, and the expression territory
    (cell) from the TSV rows. Columns with NaN values are dropped.

    Parameters:
    - filepath (str): Path to the TSV file.
    - mapping_filepath (str): Path to the TSV file containing the mapping between KH2012 and KY21.

    Returns:
    - pd.DataFrame: A preprocessed DataFrame.
    """

    # Read the TSV file into a DataFrame
    df = pd.read_csv(
        filepath, sep="\t", header=None, names=["URL", "Stage", "Gene", "Territory"]
    )

    # Extract stage in parentheses
    df["Stage"] = df["Stage"].str.extract(r"\((.*?)\)", expand=False)

    # Extract the KH number using regex
    df["Gene"] = df["Gene"].str.extract(r"(KH2012:KH\.[A-Z]\d+\.\d+)", expand=False)

    # Extract the expression territory using regex
    df["Territory"] = df["Territory"].str.extract(r"([ABab]\d+\.\d+)", expand=False)

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


def pad_compatible(mat: pd.DataFrame, adata: anndata.AnnData) -> pd.DataFrame:
    """
    Adjust the adjacency matrix to have the same number and order of columns as the AnnData object.

    Parameters:
    - mat (pd.DataFrame): Adjacency matrix with cells in rows and genes in columns.
    - adata (anndata.AnnData): An AnnData object.

    Returns:
    - pd.DataFrame: Adjusted adjacency matrix.
    """

    # Extract the gene names from the AnnData object's variable names
    genes = adata.var_names.tolist()

    # Reindex the adjacency matrix to match the AnnData object's genes
    # This will introduce NaN for missing genes which should be filled with zeros
    mat_padded = mat.reindex(columns=genes).fillna(0).astype(int)

    return mat_padded


def normalize_cols(mat: pd.DataFrame) -> pd.DataFrame:
    """
    Normalize all column vectors in a matrix using L2 norm.

    Parameters:
    - mat (pd.DataFrame): Input matrix.

    Returns:
    - pd.DataFrame: Matrix with normalized column vectors.
    """

    # Calculate the L2 norm for each column
    l2_norms = np.linalg.norm(mat.values, axis=0)

    # Normalize each column by its L2 norm
    normalized_mat = mat.divide(l2_norms, axis=1)

    return normalized_mat


def process_and_normalize(
    df: pd.DataFrame, adata: anndata.AnnData
) -> Dict[str, pd.DataFrame]:
    """
    Split the DataFrame by stage, pad each split with the AnnData object,
    and return a dictionary of normalized matrices.

    Parameters:
    - df (pd.DataFrame): DataFrame returned by preprocess_tsv.
    - adata (anndata.AnnData): An AnnData object.

    Returns:
    - Dict[str, pd.DataFrame]: Dictionary with stage as keys and normalized matrices as values.
    """

    result_dict = {}

    # Split DataFrame based on 'Stage' and process each split
    for stage, stage_df in df.groupby("Stage"):
        # Pad the DataFrame split with the AnnData object
        padded_matrix = pad_compatible(build_from_df(stage_df), adata)

        # Normalize the padded matrix
        normalized_matrix = normalize_cols(padded_matrix)

        # Store in the result dictionary
        result_dict[stage] = normalized_matrix

    return result_dict
