# ==============================================================================
# Project: Preterm EV Multi-Omics - Proteomics
# Script:  03_prot_enrich_compartments.R
# Author:  Nicole Grinovero
# Date:    Jan 2026
#
# Aim:
#   Perform subcellular localization enrichment (Jensen_COMPARTMENTS) on the
#   set of measured proteins (gene symbols). Saves the top N terms shown in
#   Fig. 3b and the full list of terms passing the filters (Supplementary Table S1).
#
# Input:
#   - data/proteomics/processed/proteins_norm_imputed.txt
#
# Output:
#   - outputs/tables/preterm/proteomics/enrich_Jensen_COMPARTMENTS.txt      (top N, Fig. 3b)
#   - outputs/tables/preterm/proteomics/enrich_Jensen_COMPARTMENTS_all.txt  (Suppl. Table S1)
#
# Notes:
#   - This is a descriptive enrichment on the measured set (not differential).
#   - Enrichr requires gene symbols; we clean and de-duplicate them.
#   - Enrichr tests the list against its own library background; results depend
#     on the library version available online at the time of the query.
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(enrichR)
})

# ------------------------------------------------------------------------------
# 0) Parameters
# ------------------------------------------------------------------------------
DB_MAIN <- "Jensen_COMPARTMENTS"
FDR_CUT <- 0.01
MIN_OR <- 1
TOP_N <- 10
S1_MAX_TERMS <- 50    # terms reported in Suppl. Table S1 (as published; Inf = all terms passing the filters)

# ------------------------------------------------------------------------------
# 1) Paths
# ------------------------------------------------------------------------------
prot_file <- here("data", "proteomics", "processed", "proteins_norm_imputed.txt")

out_dir <- here("outputs", "tables", "preterm", "proteomics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_file <- file.path(out_dir, "enrich_Jensen_COMPARTMENTS.txt")
out_all  <- file.path(out_dir, "enrich_Jensen_COMPARTMENTS_all.txt")

# ------------------------------------------------------------------------------
# 2) Load data and extract gene symbols
# ------------------------------------------------------------------------------
expr <- read_delim(prot_file, delim = "\t", show_col_types = FALSE)
stopifnot("GeneName" %in% names(expr))

symbols <- expr %>%
  filter(!is.na(GeneName)) %>%
  distinct(GeneName) %>%
  pull(GeneName)

cat("Gene symbols used for Enrichr:", length(symbols), "\n")

# ------------------------------------------------------------------------------
# 3) Run Enrichr
# ------------------------------------------------------------------------------
setEnrichrSite("Enrichr")

# Check that the library exists in Enrichr (skipped if the list cannot be retrieved)
dbs <- tryCatch(listEnrichrDbs(), error = function(e) NULL)
if (!is.null(dbs) && !"libraryName" %in% names(dbs)) dbs <- NULL
if (!is.null(dbs) && !(DB_MAIN %in% dbs$libraryName)) {
  stop("Enrichr DB not found: ", DB_MAIN)
}

enr_list <- enrichr(symbols, DB_MAIN)
if (is.null(enr_list) || !(DB_MAIN %in% names(enr_list))) {
  stop("Enrichr returned no result list for DB: ", DB_MAIN)
}

enr <- enr_list[[DB_MAIN]]
if (is.null(enr) || nrow(enr) == 0) {
  stop("Enrichr returned 0 rows for DB: ", DB_MAIN)
}

# ------------------------------------------------------------------------------
# 4) Clean, filter and select TOP N terms
# ------------------------------------------------------------------------------
# Enrichr 'Overlap' is "k/K" (k = hits in the query list; K = size of the term set)
df_top <- enr %>%
  filter(Adjusted.P.value <= FDR_CUT, Odds.Ratio >= MIN_OR) %>%
  tidyr::separate(Overlap, into = c("Count", "Background"), sep = "/", convert = TRUE) %>%
  mutate(
    Term = str_replace_all(Term, "_", " "),
    Term = str_squish(Term)
  ) %>%
  arrange(Adjusted.P.value, desc(Count)) %>%
  slice_head(n = TOP_N) %>%
  select(
    Term,
    Count,
    Background,
    Odds.Ratio,
    Adjusted.P.value,
    Genes
  )


# Full list for Supplementary Table S1 (same filters and ordering as df_top)
df_all <- enr %>%
  filter(Adjusted.P.value <= FDR_CUT, Odds.Ratio >= MIN_OR) %>%
  tidyr::separate(Overlap, into = c("Count", "Background"), sep = "/", convert = TRUE) %>%
  mutate(Term = str_squish(str_replace_all(Term, "_", " "))) %>%
  arrange(Adjusted.P.value, desc(Count)) %>%   # stable: ties keep the Enrichr order
  slice_head(n = S1_MAX_TERMS) %>%
  mutate(Rank = row_number(), In_Fig3b = Rank <= TOP_N) %>%
  select(Rank, Term, Count, Background, Odds.Ratio, P.value,
         Adjusted.P.value, Combined.Score, Genes, In_Fig3b)

write_delim(df_all, out_all, delim = "\t")
cat("Terms passing filters (Suppl. Table S1):", nrow(df_all), "\n")

cat("Terms passing filters and saved:", nrow(df_top), "\n")


# ------------------------------------------------------------------------------
# 5) Save output table
# ------------------------------------------------------------------------------
write_delim(df_top, out_file, delim = "\t")
cat("Saved:", out_file, "\n")