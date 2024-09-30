# %%
import re

import numpy as np
import pandas as pd
import geopandas as gpd
from sklearn.metrics import pairwise_distances
from scipy.spatial.distance import jaccard
import seaborn as sns
import matplotlib.pyplot as plt

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
    gdfs[stage]["name"] = gdfs[stage]["name"].str.replace("*", "", regex=False)

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
        count_dfs[stage], left_on="name", right_on="Territory_eq", how="left"
    )
    for stage, _, _ in STAGES
}

# %%

fig, axs = plt.subplots(
    nrows=1, ncols=len(STAGES), figsize=(8, 4), sharey="all", dpi=300
)

for a, (k, _, l) in zip(axs, STAGES):
    v = merged_dfs[k]
    v.plot(
        column="n",
        ax=a,
        linewidth=0.8,
        categorical=False,
        vmin=1,
        vmax=20,
        missing_kwds={"color": "lightgrey"},
        cmap=sns.color_palette("rocket", as_cmap=True),
    )
    a.axis("off")
    a.set_title(k)

patch_col = axs[0].collections[0]

fig.tight_layout()
fig.colorbar(patch_col, ax=axs, shrink=0.5)
fig.savefig("marker_map.png")

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
        pairwise_distances(pattern_df, metric=jaccard),
        index=pattern_df.index,
        columns=pattern_df.index,
    )
    for stage, pattern_df in pattern_dfs.items()
}

# %%

WIDTH = [3, 3, 2, 2, 2]

for (stage, _, _), w in zip(STAGES, WIDTH):
    plt.figure(figsize=(w, 2), dpi=300)
    pattern_df = pattern_dfs[stage]
    print(np.linalg.matrix_rank(pattern_df))
    print(pattern_df.shape)

    plt.spy(pattern_df)
    plt.xlabel("Genes")
    plt.ylabel("Blastomeres")
    plt.xticks([])
    plt.yticks([])
    plt.title(stage)
    plt.tight_layout()
    plt.savefig(f"marker_{stage.replace(" ", "_")}.png")

# %%
