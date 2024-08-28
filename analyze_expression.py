"""
This module contains functions for analyzing gene expression data.
"""

from typing import Literal

import anndata
import scanpy as sc
import pandas as pd

Method = Literal["logreg", "t-test", "wilcoxon", "t-test_overestim_var"]
PostHoc = Literal["benjamini-hochberg", "bonferroni"]


def diff_expression(
    adata: anndata.AnnData,
    groupby: str = "leiden",
    method: Method = "t-test",
    posthoc: PostHoc = "benjamini-hochberg",
    top: int = 50,
    alpha: float = 0.05,
) -> pd.DataFrame:
    """
    Find differentially expressed genes for each cluster in an AnnData object.

    :param adata: The AnnData object containing the expression data.
    :type adata: anndata.AnnData
    :param groupby: The column name in `adata.obs` that contains the cluster
    labels. Default is "leiden".
    :type groupby: str, optional
    :param method: The statistical method to use for differential expression
    analysis. Default is "t-test".
    :type method: Method, optional
    :param posthoc: The post-hoc correction method to use. Default is
    "benjamini-hochberg".
    :type posthoc: PostHoc, optional
    :param top: The number of top differentially expressed genes to select for
    each cluster. Default is 50.
    :type top: int, optional
    :param alpha: The significance level for determining differential
    expression. Default is 0.05.
    :type alpha: float, optional

    :return: A dataframe containing the differentially expressed genes for each
    cluster, sorted by log-fold change.
    :rtype: pd.DataFrame
    """

    sc.tl.rank_genes_groups(adata, groupby, method=method, corr_method=posthoc)
    result = adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names

    df_diff = pd.concat(
        [
            pd.DataFrame(
                {
                    "group": group,
                    "gene": result["names"][group],
                    "p_adj": result["pvals_adj"][group],
                    "log2fc": result["logfoldchanges"][group],
                }
            )
            for group in groups
        ]
    )

    return (
        df_diff.groupby("group")
        .apply(
            lambda g: g[g["p_adj"] < alpha]
            .sort_values(by="log2fc", ascending=False)
            .head(top)
        )
        .reset_index(drop=True)
    )
