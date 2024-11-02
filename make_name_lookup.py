import pandas as pd

df = pd.read_csv("test2.gtf", sep="\t", header=None, comment="#")

df["gene_id"] = df[8].str.extract(r'gene_id "([^"]+)"')
df["gene_name"] = df[8].str.extract(r'gene_name "([^"]+)"')

df[["gene_id", "gene_name"]].dropna().drop_duplicates().to_csv(
    "ky2021_gene_names.tsv", sep="\t", index=False
)
