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
from liana.mt import rank_aggregate

# figure settings
sc.settings.set_figure_params(dpi=200, frameon=False)
sc.set_figure_params(dpi=200, facecolor="white")
sc.set_figure_params(figsize=(5, 5))

import os
print(os.getcwd())
print(sorted(os.listdir()))

import glob
print(glob.glob("*.h5ad"))

# https://www.sc-best-practices.org/mechanisms/cell_cell_communication.html

import liana as li
import scanpy as sc

# RESOURCES
import liana as li
import pandas as pd

resources = li.resource.show_resources()

summary = []
for r in resources:
    res = li.resource.select_resource(r)
    n_total  = res.shape[0]
    n_unique = res[['ligand', 'receptor']].drop_duplicates().shape[0]
    summary.append({
        'resource'            : r,
        'total_rows'          : n_total,
        'unique_interactions' : n_unique
    })

df_summary = pd.DataFrame(summary).sort_values('unique_interactions', ascending=False)
print(df_summary)

res = li.resource.select_resource('cellphonedb')
print(res.head())

# METHODS
li.mt.show_methods()

# LIANA will by default use the .raw attribute of AnnData.
# If you wish to use .X set use_raw to False, or specify a layer.
# LIANA will also by default use the 'consensus' resource to infer ligand-receptor interactions.
# This resource was created as a consensus from the resources literature-curated resources in OmniPath,
# and uses human gene symbols.

adata = sc.datasets.pbmc68k_reduced()
print(adata)
print(adata.obs.head(3))
print(adata.obsm)
print(adata.obsm["X_pca"])
print(adata.obsm["X_umap"])
print(adata.var.head(3))
print(adata.varm)
print(adata.varm["PCs"])

sc.pl.umap(adata, color='bulk_labels', title='', frameon=False)
plt.show()

# check if raw exists first
print(adata.raw)
print(adata.raw.X[1:3, 1:3])
print(adata.raw.X[1:3, 1:3].toarray())
print(adata.raw.var_names[:10].tolist())
print("Raw shape:", adata.raw.X.shape)

rank_aggregate.describe()

print('''
| Level                      | What happens                                                        |
| -------------------------- | ------------------------------------------------------------------- |
| **Resource (`consensus`)** | curated integration of LR pairs across databases (NOT simple union) |
| **Scoring (RRA)**          | statistical consensus across methods                                |
''')

print("Using Consensus db")

res = li.resource.select_resource('consensus')
print(res.head())
print(res[['ligand', 'receptor']].drop_duplicates().shape[0])

res = li.resource.select_resource('consensus')

# import all individual methods
from liana.method import singlecellsignalr, connectome, cellphonedb, natmi, logfc, cellchat, geometric_mean

# run cellphonedb
cellphonedb(adata,
            groupby       = 'bulk_labels',
            resource_name = 'consensus',
            expr_prop     = 0.1,
            verbose       = True,
            key_added     = 'cpdb_res')

# by default, liana's output is saved in place:
print(adata.uns['cpdb_res'].head())

print('''
LIANA / CellPhoneDB Metrics Overview

Expression filtering
*_props
Proportion of cells expressing a given gene within a cell type.
By default, interactions are filtered out if either ligand or receptor is expressed in <10% of cells,
based on the assumption that meaningful cell-cell communication requires sufficient expression within a population.

Expression magnitude
*_means : Mean expression of a gene within a cell type.
* lr_means : Mean expression of the ligand-receptor pair.
Used as a proxy for interaction strength (magnitude) - higher values indicate stronger interactions.

Interaction specificity
* cellphone_pvals
Permutation-based p-values estimating how specific an interaction is between two cell types.
Lower values indicate higher specificity (i.e., less likely to occur by chance).
''')

print(list(adata.uns.keys()))

print("VISUALIZATION:")

li.pl.dotplot(adata         = adata,
              colour        = 'lr_means',
              size          = 'cellphone_pvals',
              inverse_size  = True,
              source_labels = ['CD34+', 'CD56+ NK', 'CD14+ Monocyte'],
              target_labels = ['CD34+', 'CD56+ NK'],
              figure_size   = (8, 7),
              filter_fun    = lambda x: x['cellphone_pvals'] <= 0.05,
              uns_key       = 'cpdb_res')
plt.show()

print("RANK AGGREGATE:")

rank_aggregate.describe()

print('''
| rank_aggregate | Meaning                           |
| -------------- | --------------------------------- |
| ~0 - 0.01      | Very strong / highly supported    |
| ~0.01 - 0.05   | Strong                            |
| ~0.05 - 0.2    | Moderate                          |
| >0.2           | Weak / low confidence             |
''')

print('''
| magnitude_rank | specificity_rank | Interpretation                |
| -------------- | ---------------- | ----------------------------- |
| low            | low              | strong AND specific (best)    |
| low            | high             | strong but not specific       |
| high           | low              | specific but weak             |
| high           | high             | weak + nonspecific (ignore)   |
''')

# Run rank_aggregate
li.mt.rank_aggregate(adata,
                     groupby       = 'bulk_labels',
                     resource_name = 'consensus',
                     expr_prop     = 0.1,
                     verbose       = True)

print(adata.uns['liana_res'].head())

print(list(adata.uns.keys()))

liana_res = adata.uns["liana_res"]
print(liana_res.head())
print(adata.uns["bulk_labels_colors"])
print(adata.uns["neighbors"])

for key in adata.uns.keys():
    val = adata.uns[key]
    print(f"\n{key}: {type(val).__name__}", end=" ")
    if hasattr(val, "shape"):
        print(f"shape={val.shape}")
    elif isinstance(val, dict):
        print(f"keys={list(val.keys())}")
    elif isinstance(val, list):
        print(f"len={len(val)}")
    else:
        print(f"= {val}")

# Best interactions
best = liana_res[
    (liana_res["magnitude_rank"]   <= 0.05) &
    (liana_res["specificity_rank"] <= 0.05)
].sort_values("magnitude_rank")

print(f"Top hits: {len(best)}")
print(best[["source", "target", "ligand_complex", "receptor_complex",
            "magnitude_rank", "specificity_rank"]].head(20))

import warnings
pd.options.mode.chained_assignment = None

li.pl.dotplot(
    adata             = adata,
    colour            = "magnitude_rank",
    size              = "specificity_rank",
    inverse_size      = True,
    inverse_colour    = True,
    source_labels     = ["CD34+", "CD56+ NK", "CD14+ Monocyte"],
    target_labels     = ["CD34+", "CD56+ NK"],
    top_n             = 10,
    orderby           = "magnitude_rank",
    orderby_ascending = True,
    figure_size       = (8, 7)
)
plt.show()

my_plot = li.pl.dotplot(adata          = adata,
                        colour         = 'magnitude_rank',
                        inverse_colour = True,
                        size           = 'specificity_rank',
                        inverse_size   = True,
                        source_labels  = ['CD34+', 'CD56+ NK', 'CD14+ Monocyte'],
                        target_labels  = ['CD34+', 'CD56+ NK'],
                        filter_fun     = lambda x: x['specificity_rank'] <= 0.01)
print(my_plot)

print("Circle Plot:")

li.pl.circle_plot(adata,
                  groupby       = 'bulk_labels',
                  score_key     = 'magnitude_rank',
                  inverse_score = True,
                  source_labels = 'CD34+',
                  filter_fun    = lambda x: x['specificity_rank'] <= 0.05,
                  pivot_mode    = 'counts',
                  figure_size   = (10, 10))
plt.show()

print("Customizing LIANA's RANK AGGREGATE: the user can choose to include only a subset of the methods")

methods = [logfc, geometric_mean]
new_rank_aggregate = li.mt.AggregateClass(li.mt.aggregate_meta, methods=methods)

new_rank_aggregate(adata,
                   groupby   = 'bulk_labels',
                   expr_prop = 0.1,
                   verbose   = True,
                   n_perms   = None,
                   use_raw   = True,
                   key_added = "rank_agg_res")

print(adata.uns['rank_agg_res'].head())

from liana.method import rank_aggregate

rank_aggregate(
    adata,
    groupby   = "bulk_labels",
    expr_prop = 0.1,
    verbose   = True,
    n_perms   = None,
    use_raw   = True,
    key_added = "rank_agg_res"
)

print(adata.uns["rank_agg_res"].head())

print('''
What is the TRI-MEAN ?

The trimean is a robust measure of central tendency that combines:
* median
* lower quartile (Q1)
* upper quartile (Q3)

import numpy as np

def trimean(x):
    q1     = np.percentile(x, 25)
    median = np.percentile(x, 50)
    q3     = np.percentile(x, 75)
    return (q1 + 2 * median + q3) / 4

Interpretation:
Gives 50% weight to the median
Gives 25% weight to Q1
Gives 25% weight to Q3
''')

print('''
Trimean = a robust average combining median and quartiles, used to stabilize gene expression estimates in noisy single-cell data.

| Metric      | Behavior                       |
| ----------- | ------------------------------ |
| Mean        | sensitive to outliers          |
| Median      | ignores distribution shape     |
| Trimean     | robust + still reflects spread |
''')

print('''
CellChat computes interaction probability as:

def cellchat_prob(L_star, R_star, Kh=0.5):
    numerator   = L_star * R_star
    denominator = Kh + L_star * R_star
    return numerator / denominator

Where:
L* = trimean expression of the ligand in the source cell type
R* = trimean expression of the receptor in the target cell type
Kh = Hill normalization constant (default = 0.5)

This is a Hill function - it saturates at 1 when expression is very high, and approaches 0 when expression is low.
The Kh constant controls the steepness of the curve.
''')

li.rs.describe_metalinks()

print('''
Annotating Ligand-Receptors:

We use commonly PROGENy pathway weights;
Disease Annotations: DisGeNet;
Enrichment score: Decoupler: https://decoupler.readthedocs.io/en/latest/
Intracellular Signaling: a protein-protein interaction network;
Transcription Factor Regulons;
Metabolite-Receptor Interaction.
''')

res = li.resource.get_metalinks(
    source           = ['Stich', 'CellPhoneDB', 'NeuronChat'],
    tissue_location  = 'Brain'
)

res_neuronchat = res[res['source'] == 'NeuparonChat']
print(res_neuronchat.head())
print(res_neuronchat.shape)
print(res_neuronchat.tail())