# %%
import re

import numpy as np
import pandas as pd
import geopandas as gpd
from sklearn.metrics import pairwise_distances
from scipy.spatial.distance import jaccard
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects

# %%
PATTERN_STAGES = r"(early|mid|late) (gastrula|neurula)"
PATTERN_CELLS = r"[Aa]\d+\.\d+$"

STAGES = [
    ("mid gastrula", "mid_gastrula.geojson", 6),
    ("late gastrula", "late_gastrula.geojson", 7),
    ("early neurula", "early_neurula.geojson", 10),
    ("mid neurula", "mid_neurula.geojson", 11),
    ("late neurula", "late_neurula.geojson", 12),
]

# %%
df = pd.read_csv("pass_02.tsv", sep="\t")

gdfs = {stage: gpd.read_file(file) for stage, file, _ in STAGES}

for stage in gdfs:
    gdfs[stage]["name_tmp"] = gdfs[stage]["name"].str.replace("*", "", regex=False)

# %%
df["Stage"] = df["Stage"].apply(lambda x: (" ".join(re.findall(PATTERN_STAGES, x)[0])))
df = df[df["Territory_eq"].str.contains(PATTERN_CELLS)]
df = df[["Stage", "Gene", "Territory_eq"]].drop_duplicates()

dfs = dict(tuple(df.groupby("Stage")))

count_dfs = {
    stage: df.groupby("Territory_eq")["Gene"]
    .nunique()
    .reset_index()
    .rename(columns={"Gene": "n"})
    for stage, df in dfs.items()
}

# %%
merged_dfs = {
    stage: gdfs[stage].merge(
        count_dfs[stage], left_on="name_tmp", right_on="Territory_eq", how="left"
    )
    for stage, _, _ in STAGES
}

# %%

for stage, _, _ in STAGES:
    fig, ax = plt.subplots(figsize=(6.5, 6.5), dpi=300)
    v = merged_dfs[stage]
    v.plot(
        column="n",
        ax=ax,
        linewidth=0.8,
        categorical=False,
        # vmin=1,
        # vmax=20,
        missing_kwds={"color": "lightgrey"},
        cmap=sns.color_palette("rocket", as_cmap=True),
    )
    # Add blastomere names at centroid positions
    for idx, row in v.iterrows():
        if row["geometry"] is not None and hasattr(row["geometry"], "centroid"):
            centroid = row["geometry"].centroid
            txt = ax.annotate(
                row["name"],
                (centroid.x, centroid.y),
                color="white",
                fontsize=6,
                ha="center",
                va="center",
            )
            txt.set_path_effects(
                [
                    path_effects.Stroke(linewidth=1.5, foreground="black"),
                    path_effects.Normal(),
                ]
            )
    ax.axis("off")
    ax.set_title(stage, fontsize=10)

    patch_col = ax.collections[0]
    fig.colorbar(patch_col, ax=ax, shrink=0.5)
    fig.savefig(f"marker_map_{stage.replace(' ', '_')}.png")


# %%
pattern_dfs = {
    stage: pd.pivot_table(
        df,
        values="Gene",
        index="Territory_eq",
        columns="Gene",
        aggfunc="size",
        fill_value=0,
    ).astype(bool)
    for stage, df in dfs.items()
}

# %%
jaccard_dfs = {
    stage: pd.DataFrame(
        1 - pairwise_distances(pattern_df, metric=jaccard),
        index=pattern_df.index,
        columns=pattern_df.index,
    )
    for stage, pattern_df in pattern_dfs.items()
}


# %%

for stage, _, _ in STAGES:
    plt.figure(figsize=(6.5, 6.5), dpi=300)
    pattern_df = pattern_dfs[stage]
    print(np.linalg.matrix_rank(pattern_df))
    print(pattern_df.shape)

    plt.spy(pattern_df)
    plt.xlabel("Genes")
    plt.ylabel("Blastomeres")
    plt.title(stage)
    plt.tight_layout()

    plt.yticks(
        ticks=np.arange(len(pattern_df.index)), labels=pattern_df.index, fontsize=5
    )
    plt.xticks(
        ticks=np.arange(len(pattern_df.columns)),
        labels=pattern_df.columns,
        fontsize=5,
        rotation=90,
    )

    plt.savefig(f"marker_mat_{stage.replace(' ', '_')}.png")

# %%
for stage, _, _ in STAGES:
    jaccard_df = jaccard_dfs[stage]
    g = sns.clustermap(
        jaccard_df,
        cmap="rocket",
        figsize=(6.5, 6.5),
        row_cluster=True,
        col_cluster=True,
        cbar_kws={"shrink": 0.5},
        xticklabels=True,
        yticklabels=True,
    )
    g.ax_heatmap.set_title(f"Jaccard Similarity - {stage}")
    g.ax_heatmap.set_xlabel("Blastomeres")
    g.ax_heatmap.set_ylabel("Blastomeres")
    g.ax_heatmap.set_xticklabels(
        g.ax_heatmap.get_xticklabels(),
        fontsize=5,
        rotation=90,
    )
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=5)
    plt.tight_layout()
    plt.savefig(
        f"jaccard_{stage.replace(' ', '_')}.png",
        dpi=300,
    )

# %%
