# ============================================================
# Extract genes from CellChatDB interaction names
# Validate gene symbols (HGNC / human)
# Save validated genes + suggested fixes (CSV + TXT)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(HGNChelper)
})

# ---- Input ----
# interactions <- CellChatDB$interaction
interactions = read.delim("CellChatDB_interaction.txt", stringsAsFactors = FALSE, sep="\t", header=TRUE)

stopifnot("interaction_name" %in% colnames(interactions))

# ---- 1) Extract genes from interaction_name ----
genes_long <- interactions %>%
  select(interaction_name) %>%
  mutate(gene = strsplit(as.character(interaction_name), "_")) %>%
  unnest(gene) %>%
  mutate(gene = trimws(gene)) %>%
  filter(!is.na(gene), gene != "")

unique_genes <- sort(unique(genes_long$gene))

cat("Total gene tokens:", nrow(genes_long), "\n")
cat("Unique genes:", length(unique_genes), "\n")

# ---- 2) Validate gene symbols (HGNC) ----
gene_check <- checkGeneSymbols(unique_genes)

cat("\nValidation summary:\n")
print(
  gene_check %>%
    count(Approved) %>%
    mutate(percent = round(100 * n / sum(n), 2))
)

# ---- 3) Extract suggested fixes for invalid symbols ----
suggested_fixes <- gene_check %>%
  filter(!Approved, !is.na(Suggested.Symbol)) %>%
  distinct(x, Suggested.Symbol) %>%
  arrange(x)

cat("\nInvalid symbols with suggested fixes (first 30):\n")
print(head(suggested_fixes, 30))

# ---- 4) Save outputs ----
outdir <- "CellChatDB_gene_validation"
dir.create(outdir, showWarnings = FALSE)

# Full validation table
write.csv(
  gene_check,
  file = file.path(outdir, "interaction_genes_validated.csv"),
  row.names = FALSE
)

# Suggested fixes (CSV)
write.csv(
  suggested_fixes,
  file = file.path(outdir, "interaction_genes_suggested_fixes.csv"),
  row.names = FALSE
)

# Suggested fixes (TXT)
write.table(
  suggested_fixes,
  file = file.path(outdir, "interaction_genes_suggested_fixes.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Valid genes (TXT)
write.table(
  gene_check$x[gene_check$Approved],
  file = file.path(outdir, "interaction_genes_valid.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# Invalid genes (TXT)
write.table(
  gene_check$x[!gene_check$Approved],
  file = file.path(outdir, "interaction_genes_invalid.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

cat("\nSaved outputs to:", normalizePath(outdir), "\n")
cat(" - interaction_genes_validated.csv\n")
cat(" - interaction_genes_suggested_fixes.csv\n")
cat(" - interaction_genes_suggested_fixes.txt\n")
cat(" - interaction_genes_valid.txt\n")
cat(" - interaction_genes_invalid.txt\n")
