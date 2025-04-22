import re

import pandas as pd

# Generate the BLAST map from the following:
# \time blastp -db swissprot -taxids 9606 -query ky2021p.fasta -parse_deflines \
# -outfmt "7 qacc sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
# -out ky2021_swissprot.txt -num_threads 64 -mt_mode 1

# \time blastp -db ky2021p.fasta -taxids 1774208 -query kh2012p.fasta \
# -parse_deflines \
# -outfmt "7 qacc sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" -out kh2012_ky2021.txt \
# -num_threads 64 -mt_mode 1

# \time blastp -db ky2021p.fasta -taxids 1774208 -query kh2013p.fasta \
# -parse_deflines \
# -outfmt "7 qacc sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" -out kh2013_ky2021.txt \
# -num_threads 64 -mt_mode 1

df_ky_sp = pd.read_csv(
    "ky2021_swissprot.txt",
    sep="\t",
    names=[
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "coverage",
    ],
    comment="#",
)
df_ky_sp["qseqid"] = df_ky_sp["qseqid"].apply(lambda x: re.sub(r"\.v.+", "", x))
df_ky_sp["sseqid"] = (
    df_ky_sp["sseqid"]
    .apply(lambda x: re.sub(r"sp\|", "", x))
    .apply(lambda x: re.sub(r"\..+", "", x))
)

df_ky_sp = df_ky_sp[df_ky_sp["evalue"] < 0.05]

df_ky_sp = df_ky_sp.loc[df_ky_sp.groupby("qseqid")["evalue"].idxmin()]

df_ky_sp = df_ky_sp.rename(
    columns={
        "sseqid": "uniprot",
        "qseqid": "KY2021",
    },
)

df_kh2012 = pd.read_csv(
    "kh2012_ky2021.txt",
    sep="\t",
    names=[
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "coverage",
    ],
    comment="#",
)

df_kh2012["qseqid"] = df_kh2012["qseqid"].apply(lambda x: re.sub(r"\.v.+", "", x))
df_kh2012["sseqid"] = df_kh2012["sseqid"].apply(lambda x: re.sub(r"\.v.+", "", x))

df_kh2012 = df_kh2012[df_kh2012["evalue"] < 0.05]

df_kh2012 = df_kh2012.loc[df_kh2012.groupby("qseqid")["evalue"].idxmin()]

df_kh2012 = df_kh2012.rename(
    columns={
        "qseqid": "KH2012",
        "sseqid": "KY2021",
    },
)[["KH2012", "KY2021"]]

df_kh2012.to_csv(
    "kh2012_ky2021_map.tsv",
    sep="\t",
    index=False,
)

df_kh2013 = pd.read_csv(
    "kh2013_ky2021.txt",
    sep="\t",
    names=[
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "coverage",
    ],
    comment="#",
)

df_kh2013["qseqid"] = df_kh2013["qseqid"].apply(lambda x: re.sub(r"\.v.+", "", x))
df_kh2013["sseqid"] = df_kh2013["sseqid"].apply(lambda x: re.sub(r"\.v.+", "", x))

df_kh2013 = df_kh2013[df_kh2013["evalue"] < 0.05]

df_kh2013 = df_kh2013.loc[df_kh2013.groupby("qseqid")["evalue"].idxmin()]

df_kh2013 = df_kh2013.rename(
    columns={
        "qseqid": "KH2013",
        "sseqid": "KY2021",
    },
)[["KH2013", "KY2021"]]

df_kh2013.to_csv(
    "kh2013_ky2021_map.tsv",
    sep="\t",
    index=False,
)

df_ky_sp = df_ky_sp.merge(
    df_kh2012,
    how="left",
    left_on="KY2021",
    right_on="KY2021",
).merge(
    df_kh2013,
    how="left",
    left_on="KY2021",
    right_on="KY2021",
)

(
    pd.merge(
        df_ky_sp,
        pd.read_csv("uniprot_data.csv"),
        left_on="uniprot",
        right_on="uniprot",
    )
    .groupby("KY2021")
    .apply(lambda x: x.nsmallest(1, "evalue"))
    .reset_index(drop=True)[
        [
            "KY2021",
            "KH2012",
            "KH2013",
            "uniprot",
            "fullname",
        ]
    ]
).to_csv(
    "ky2021_swissprot_map.csv",
    index=False,
)
