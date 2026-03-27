#!/usr/bin/env python

import decoupler as dc
import liana as li
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import session_info

from liana.method import cellphonedb

# figure settings
sc.settings.set_figure_params(dpi=200, frameon=False)
sc.set_figure_params(dpi=200, facecolor="white")
sc.set_figure_params(figsize=(5, 5))

import os
print(os.getcwd())
print(sorted(os.listdir()))

# Read in
adata = sc.read("kang.h5ad")
print(adata)

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

print(adata.X[1:5, 1:5].toarray())

# Store the counts for later use
adata.layers["counts"] = adata.X.copy()

# Rename label to condition, replicate to patient
adata.obs = adata.obs.rename({"label": "condition", "replicate": "patient"}, axis=1)

# assign sample
adata.obs["sample"] = (
    adata.obs["condition"].astype("str") + "&" + adata.obs["patient"].astype("str")
)

print(adata.obs.head(5))
print(adata.var.head(5))

print(adata.obs["condition"].unique())
print(adata.obs["cell_type"].unique())
print(adata.obs["condition"].value_counts())
print(adata.obs["cell_type"].value_counts())

# log1p normalize the data
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

sc.pl.umap(adata, color=["condition", "cell_type"], frameon=False)
plt.show()

# Unique values
for col in ["condition", "cluster", "cell_type", "patient"]:
    print(f"\n{col}: {adata.obs[col].unique().tolist()}")

print("\n===== Cells per condition =====")
print(adata.obs["condition"].value_counts())

print("\n===== Patients per condition =====")
print(adata.obs.groupby("condition")["patient"].nunique())

print("\n===== Cells per patient per condition =====")
print(adata.obs.groupby(["condition", "patient"]).size().reset_index(name="n_cells"))

# Cells per cell_type per patient
print("\n===== Cells per cell_type x patient =====")
tbl = (adata.obs
       .groupby(["cell_type", "patient"])
       .size()
       .reset_index(name="n_cells")
       .pivot(index="cell_type", columns="patient", values="n_cells")
       .fillna(0)
       .astype(int))
tbl["TOTAL"] = tbl.sum(axis=1)
tbl = tbl.sort_values("TOTAL", ascending=False)
print(tbl.to_string())

print("\n===== Cells per cell_type x condition =====")
tbl2 = (adata.obs
        .groupby(["cell_type", "condition"])
        .size()
        .reset_index(name="n_cells")
        .pivot(index="cell_type", columns="condition", values="n_cells")
        .fillna(0)
        .astype(int))
tbl2["TOTAL"] = tbl2.sum(axis=1)
tbl2 = tbl2.sort_values("TOTAL", ascending=False)
print(tbl2.to_string())

print("\n===== Cells per cell_type x patient x condition =====")
tbl3 = (adata.obs
        .groupby(["cell_type", "patient", "condition"])
        .size()
        .reset_index(name="n_cells"))
print(tbl3.to_string())

# ── LIGAND RECEPTOR inference ─────────────────────────────────
print("LIGAND RECEPTOR inference :")

adata_stim = adata[adata.obs["condition"] == "stim"].copy()
print(adata_stim)

cellphonedb(
    adata_stim,
    groupby        = "cell_type",
    use_raw        = False,
    return_all_lrs = True,
    verbose        = True
)

print(adata_stim.uns["liana_res"].head())

li.rs.show_resources()
li.method.show_methods()

# ── Visual Exploration ────────────────────────────────────────
print("Visual Exploration")

li.pl.dotplot(
    adata             = adata_stim,
    colour            = "lr_means",
    size              = "cellphone_pvals",
    inverse_size      = True,
    source_labels     = ["CD4 T cells", "B cells", "FCGR3A+ Monocytes"],
    target_labels     = ["CD8 T cells", "CD14+ Monocytes", "NK cells"],
    orderby           = "lr_means",
    orderby_ascending = False,
    top_n             = 20,
    figure_size       = (10, 14),
    size_range        = (1, 6),
)
plt.show()

# ── Generating a Ligand-Receptor consensus with LIANA ─────────
print("Generating a Ligand-Receptor consensus with LIANA")

from liana.method import rank_aggregate

rank_aggregate(
    adata_stim,
    groupby        = "cell_type",
    return_all_lrs = True,
    use_raw        = False,
    verbose        = True
)

print(adata_stim.uns["liana_res"].drop_duplicates(
    ["ligand_complex", "receptor_complex"]
).head())

# ── Filter and plot rank_aggregate ────────────────────────────
liana_res      = adata_stim.uns["liana_res"].copy()
liana_filtered = liana_res.loc[liana_res["specificity_rank"] <= 0.05].copy()
adata_stim.uns["liana_res"] = liana_filtered

li.pl.dotplot(
    adata             = adata_stim,
    colour            = "magnitude_rank",
    size              = "specificity_rank",
    inverse_colour    = True,
    inverse_size      = True,
    source_labels     = ["CD4 T cells", "B cells", "FCGR3A+ Monocytes"],
    target_labels     = ["CD8 T cells", "CD14+ Monocytes", "NK cells"],
    orderby           = "magnitude_rank",
    orderby_ascending = True,
    top_n             = 20,
    figure_size       = (10, 16),
    size_range        = (1, 6),
)
plt.show()

adata_stim.uns["liana_res"] = liana_res