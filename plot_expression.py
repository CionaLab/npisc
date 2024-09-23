# %%
import re

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

fig, axs = plt.subplots(nrows=1, ncols=len(STAGES), figsize=(8, 4), sharey="all")

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

plt.tight_layout()
plt.colorbar(patch_col, ax=axs, shrink=0.5)
plt.savefig("marker_map.png", dpi=300)

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
for stage, _, _ in STAGES:
    pattern_df = pattern_dfs[stage]
    jaccard_df = jaccard_dfs[stage]

    g = sns.clustermap(pattern_df, cmap="cool")
    g.figure.suptitle(f"Pattern Heatmap for {stage}")
    g.figure.subplots_adjust(top=0.95, right=0.8)
    g.ax_cbar.set_position((0.9, 0.2, 0.03, 0.4))

    g = sns.clustermap(jaccard_df, method="complete", metric="jaccard", cmap="cool")
    g.figure.suptitle(f"Jaccard Distance Heatmap for {stage}")
    g.figure.subplots_adjust(top=0.95, right=0.8)
    g.ax_cbar.set_position((0.9, 0.2, 0.03, 0.4))

    plt.show()

# %%
