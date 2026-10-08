# ==============================================================================
# Project: Preterm EV Multi-Omics – Proteomics
# Script:  01_prot_cleaning.R
# Author:  Nicole Grinovero
# Date:    Jan 2026
#
# Aim:
#   Clean the Spectronaut protein pivot export and generate a curated protein
#   matrix (annotations + samples) ready for normalization/imputation.
#
# Inputs:
#   - data/metadata/sample_metadata.txt
#   - data/proteomics/raw/20251106_112203_Andreato_Plasma_Prematuri_MagNet_Report.tsv
#     (Spectronaut report; protein group quantities exported on log2 scale)
#
# Output:
#   - data/proteomics/processed/proteins_cleaned.txt
#
# Notes:
#   - Spectronaut often stores multiple identifiers in a single field separated
#     by ';' (e.g. multiple DBs / multiple accessions). We prioritize UniProt
#     SwissProt ("sp") IDs when present.
#   - Sample columns are renamed from technical raw IDs to biological sample IDs
#     using the metadata table.
#   - Quantities are already log2-transformed in the Spectronaut export; no log
#     transformation is applied here or downstream.
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
})

# ------------------------------------------------------------------------------
# 0) Helper: extract preferred ID from Spectronaut multi-ID fields
# ------------------------------------------------------------------------------
# Spectronaut fields can look like:
#   PG.Database:       "sp;tr;sp"
#   PG.ProteinGroups:  "P12345;Q99999;P67890"
# We keep the first ID belonging to target_db ("sp") if available, otherwise
# we fall back to the first ID in the list.
extract_id <- function(db_str, id_str, target_db = "sp") {
  if (is.na(db_str) || is.na(id_str)) return(NA_character_)
  
  db_values <- str_split(db_str, ";", simplify = TRUE) |> str_trim()
  id_values <- str_split(id_str, ";", simplify = TRUE) |> str_trim()
  
  if (length(id_values) == 1) return(id_values[1])
  
  idx <- which(db_values == target_db)
  if (length(idx) > 0) return(id_values[idx[1]])
  
  id_values[1]
}

# ------------------------------------------------------------------------------
# 1) Paths
# ------------------------------------------------------------------------------
path_meta <- here("data", "metadata", "sample_metadata.txt")
path_raw  <- here(
  "data", "proteomics", "raw",
  "20251106_112203_Andreato_Plasma_Prematuri_MagNet_Report.tsv"
)

out_dir  <- here("data", "proteomics", "processed")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, "proteins_cleaned.txt")

# ------------------------------------------------------------------------------
# 2) Load metadata + raw Spectronaut export
# ------------------------------------------------------------------------------
meta_prot <- read_delim(path_meta, delim = "\t", show_col_types = FALSE)
prot_raw  <- read_delim(path_raw,  delim = "\t", show_col_types = FALSE)

cat("Metadata samples:", nrow(meta_prot), "\n")
cat("Raw proteins:", nrow(prot_raw), "\n")

stopifnot(all(c("sample_id", "proteomics_raw_id") %in% names(meta_prot)))
stopifnot(!anyDuplicated(meta_prot$sample_id))

# ------------------------------------------------------------------------------
# 3) Standardize quantitative column names (technical IDs)
# ------------------------------------------------------------------------------
# Spectronaut column names can contain extra text; we extract a stable raw ID
# pattern (LPD-###-T#_#) to match the metadata proteomics_raw_id.
extract_raw_id_from_col <- function(x) {
  m <- str_match(x, "(LPD-\\d{3}-T\\d+_\\d+)")
  ifelse(!is.na(m[, 2]), m[, 2], x)
}

prot_std <- prot_raw %>%
  rename_with(.cols = where(is.numeric), .fn = extract_raw_id_from_col)

sample_cols_all <- names(prot_std) %>% str_subset("^LPD-\\d{3}-T\\d+_\\d+$")
sample_cols_raw <- intersect(sample_cols_all, meta_prot$proteomics_raw_id)

cat("Matched raw sample columns:", length(sample_cols_raw), "\n")
if (length(sample_cols_raw) == 0) {
  stop("No quantitative columns matched metadata proteomics_raw_id. Check naming pattern / metadata.")
}

# ------------------------------------------------------------------------------
# 4) Add curated identifiers (UniProtID, GeneName)
# ------------------------------------------------------------------------------
prot_clean <- prot_std %>%
  mutate(
    UniProtID = if (all(c("PG.Database", "PG.ProteinGroups") %in% names(.))) {
      map2_chr(.data$PG.Database, .data$PG.ProteinGroups, extract_id)
    } else NA_character_,
    GeneName = if (all(c("PG.Database", "PG.Genes") %in% names(.))) {
      map2_chr(.data$PG.Database, .data$PG.Genes, extract_id)
    } else NA_character_
  )

# ------------------------------------------------------------------------------
# 5) Rename sample columns (technical -> biological) and subset output columns
# ------------------------------------------------------------------------------
rename_map <- setNames(meta_prot$proteomics_raw_id, meta_prot$sample_id)

prot_renamed <- prot_clean %>%
  dplyr::rename(any_of(rename_map))

annot_cols <- c(
  "PG.ProteinDescriptions",
  "PG.ProteinGroups",
  "PG.Genes",
  "PG.MolecularFunction",
  "PG.BiologicalProcess",
  "PG.CellularComponent"
)

bio_cols <- meta_prot$sample_id

prot_final <- prot_renamed %>%
  select(
    any_of(c("UniProtID", "GeneName")),
    any_of(intersect(annot_cols, names(.))),
    all_of(intersect(bio_cols, names(.)))
  )

cat("Final table:", nrow(prot_final), "proteins x",
    sum(names(prot_final) %in% bio_cols), "samples\n")

# ------------------------------------------------------------------------------
# 6) Save
# ------------------------------------------------------------------------------
write_delim(prot_final, out_file, delim = "\t")
cat("Saved to:", out_file, "\n")
