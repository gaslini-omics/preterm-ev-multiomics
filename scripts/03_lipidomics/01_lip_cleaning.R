# ==============================================================================
# Project: Preterm EV Multi-Omics – Lipidomics
# Script:  01_lip_cleaning.R
# Author:  Nicole Grinovero
# Date:    Jan 2026
#
# Aim:
#   Clean and curate the MS-FLO-processed lipidomics area table to obtain a
#   high-quality feature-level matrix ready for normalization/imputation.
#
# Inputs:
#   - data/metadata/sample_metadata.txt
#   - data/lipidomics/raw/Area_1_2025_12_17_16_14_21.txt            (MS-DIAL export, before MS-FLO)
#   - data/lipidomics/raw/Area_1_2025_12_17_16_14_21_processed.txt  (MS-FLO processed)
#
# Output:
#   data/lipidomics/processed/lipids_cleaned.txt
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(rgoslin)
  library(LipidSigR)
})

# ------------------------------------------------------------------------------
# 0) Configuration
# ------------------------------------------------------------------------------
tp_levels <- c("T0", "T1", "T2", "T3", "T4")

THR_QC_CV <- 0.30
THR_FILL  <- 0.70

raw_dir  <- here("data", "lipidomics", "raw")
proc_dir <- here("data", "lipidomics", "processed")
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)

meta_path     <- here("data", "metadata", "sample_metadata.txt")
raw_pre_path  <- file.path(raw_dir,  "Area_1_2025_12_17_16_14_21.txt")
raw_expr_path <- file.path(raw_dir,  "Area_1_2025_12_17_16_14_21_processed.txt")
out_file      <- file.path(proc_dir, "lipids_cleaned.txt")

CONTAMINANTS_REGEX <- c(
  "^RIKEN", "from\\s+Blank", "chemical", "reagent",
  "phthalate", "siloxane", "polyethylene",
  "docosenamide", "erucamide", "oleamide"
)

info_cols <- c(
  "Alignment ID", "LipidName",
  "Average Mz", "Average Rt(min)", "Adduct type", "Ontology", "INCHIKEY", "SMILES"
)

# ------------------------------------------------------------------------------
# 0.1) Helper: recover/clean names from MS-FLO output
# ------------------------------------------------------------------------------
clean_lipid_names <- function(df, raw_pre_df) {
  stopifnot(all(c("Alignment ID", "Metabolite name", "MS/MS matched") %in% names(df)))
  stopifnot(all(c("Alignment ID", "Metabolite name") %in% names(raw_pre_df)))
  
  id2name <- raw_pre_df %>%
    transmute(
      AlignmentID   = as.character(`Alignment ID`),
      Name_Original = as.character(`Metabolite name`)
    ) %>%
    distinct(AlignmentID, .keep_all = TRUE)
  
  df %>%
    filter(str_detect(`MS/MS matched`, "True")) %>%
    mutate(
      id_chr   = as.character(`Alignment ID`),
      is_merge = str_detect(id_chr, "_"),
      id_L     = str_split_fixed(id_chr, "_", 2)[, 1],
      id_R     = str_split_fixed(id_chr, "_", 2)[, 2]
    ) %>%
    left_join(id2name, by = c("id_L" = "AlignmentID")) %>%
    left_join(id2name, by = c("id_R" = "AlignmentID"), suffix = c("_L", "_R")) %>%
    mutate(
      LipidName = case_when(
        is_merge & `MS/MS matched` == "True_False"  ~ Name_Original_L,
        is_merge & `MS/MS matched` == "False_True"  ~ Name_Original_R,
        TRUE ~ `Metabolite name`
      ),
      LipidName = str_remove(LipidName, "^.*\\|"),
      LipidName = str_remove(LipidName, "_RIKEN.*"),
      LipidName = str_remove(LipidName, "_[[:alpha:]][^_]*$"),
      LipidName = str_squish(LipidName)
    ) %>%
    select(-id_chr, -is_merge, -id_L, -id_R, -Name_Original_L, -Name_Original_R)
}

# ------------------------------------------------------------------------------
# 0.2) Helper: lipid key for safe deduplication
# ------------------------------------------------------------------------------
# Minimal rule: for lyso-lipids only, collapse "/0:0" notation
# e.g., "LPC 18:0/0:0" -> "LPC 18:0"
make_lipid_key <- function(x) {
  x <- str_squish(as.character(x))
  is_lyso <- str_detect(x, "^LP")  # LPC/LPE/LPI/LPS/...
  x[is_lyso] <- str_replace(x[is_lyso], "/0:0$", "")
  str_squish(x)
}

# ------------------------------------------------------------------------------
# 1) Load metadata + raw data
# ------------------------------------------------------------------------------
meta <- read_delim(meta_path, delim = "\t", show_col_types = FALSE)

meta_lipids <- meta %>%
  filter(!is.na(lipidomics_raw_id)) %>%
  mutate(timepoint = factor(timepoint, levels = tp_levels)) %>%
  arrange(timepoint, patient_id)

raw_pre  <- read_delim(raw_pre_path,  delim = "\t", skip = 4, show_col_types = FALSE)
raw_expr <- read_delim(raw_expr_path, delim = "\t", show_col_types = FALSE)

cat("Samples in lipidomics metadata:", nrow(meta_lipids), "\n")
cat("Raw features (input table):", nrow(raw_expr), "\n")

# ------------------------------------------------------------------------------
# 2) Lipid name curation + contaminant removal
# ------------------------------------------------------------------------------
expr_clean_names <- clean_lipid_names(raw_expr, raw_pre) %>%
  filter(!str_detect(
    LipidName,
    regex(paste(CONTAMINANTS_REGEX, collapse = "|"), ignore_case = TRUE)
  ))

cat("Features after name curation + contaminant removal:", nrow(expr_clean_names), "\n")

# ------------------------------------------------------------------------------
# 3) QC filtering (CV over QC injections)
# ------------------------------------------------------------------------------
qc_cols <- colnames(expr_clean_names) %>% str_subset("QC")
if (length(qc_cols) == 0) {
  stop("No QC columns detected (pattern 'QC'). Check input table column names.")
}

expr_qc <- expr_clean_names %>%
  rowwise() %>%
  mutate(
    qc_mean = mean(c_across(all_of(qc_cols)), na.rm = TRUE),
    qc_cv   = if_else(qc_mean > 0,
                      sd(c_across(all_of(qc_cols)), na.rm = TRUE) / qc_mean,
                      NA_real_)
  ) %>%
  ungroup() %>%
  filter(!is.na(qc_cv) & qc_cv <= THR_QC_CV)

cat("After QC filtering (QC CV <=", THR_QC_CV, "):", nrow(expr_qc), "\n")

# ------------------------------------------------------------------------------
# 4) Rename technical IDs -> biological sample IDs (metadata-driven)
# ------------------------------------------------------------------------------
rename_map <- setNames(meta_lipids$lipidomics_raw_id, meta_lipids$sample_id)
bio_cols   <- meta_lipids$sample_id

expr_renamed <- expr_qc %>%
  rename(any_of(rename_map)) %>%
  select(any_of(info_cols), all_of(bio_cols))

# Enforce sample ordering to match metadata
expr_renamed <- expr_renamed %>%
  select(any_of(info_cols), all_of(bio_cols))

# ------------------------------------------------------------------------------
# 5) Define missingness (convert non-positive values to NA)
# ------------------------------------------------------------------------------
expr_renamed <- expr_renamed %>%
  mutate(across(all_of(bio_cols), ~ if_else(. <= 0, NA_real_, as.numeric(.))))

# ------------------------------------------------------------------------------
# 6) Biological filtering (fill-rate >= THR_FILL in >=1 timepoint)
# ------------------------------------------------------------------------------
feature_stats <- expr_renamed %>%
  select(`Alignment ID`, all_of(bio_cols)) %>%
  pivot_longer(-`Alignment ID`, names_to = "sample_id", values_to = "Intensity") %>%
  left_join(meta_lipids %>% select(sample_id, timepoint), by = "sample_id") %>%
  group_by(`Alignment ID`, timepoint) %>%
  summarise(fill_rate = sum(!is.na(Intensity)) / n(), .groups = "drop_last") %>%
  summarise(max_fill = max(fill_rate), .groups = "drop")

expr_fill <- expr_renamed %>%
  left_join(feature_stats, by = "Alignment ID") %>%
  filter(max_fill >= THR_FILL)

cat("After fill-rate filtering (max fill >=", THR_FILL, "):", nrow(expr_fill), "\n")

# ------------------------------------------------------------------------------
# 7) Duplicate resolution (one feature per LipidKey)
# ------------------------------------------------------------------------------
# We first standardize a minimal key (lyso "/0:0" collapse), then keep a single
# representative feature per lipid using:
#   1) higher max_fill (less missing)
#   2) lower qc_cv   (more reproducible)
#   3) higher qc_mean (stronger signal)
expr_fill2 <- expr_fill %>%
  left_join(expr_qc %>% select(`Alignment ID`, qc_mean, qc_cv), by = "Alignment ID") %>%
  mutate(LipidKey = make_lipid_key(LipidName))

n_before <- nrow(expr_fill2)

expr_final <- expr_fill2 %>%
  group_by(LipidKey) %>%
  arrange(desc(max_fill), qc_cv, desc(qc_mean)) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(LipidName = LipidKey) %>%        # enforce standardized naming in output
  select(-LipidKey, -qc_cv, -qc_mean, -max_fill)

cat("Duplicate collapse:", n_before, "->", nrow(expr_final), "features\n")

# ------------------------------------------------------------------------------
# 8) Lipid class annotation (Goslin + LipidSigR)
# ------------------------------------------------------------------------------
parsed <- rgoslin::parseLipidNames(lipidNames = expr_final$LipidName)

goslin_status <- parsed %>%
  transmute(LipidName = Original.Name, Grammar = Grammar)

recognized <- goslin_status %>%
  filter(Grammar != "NOT_PARSEABLE") %>%
  pull(LipidName)

not_parseable <- goslin_status %>%
  filter(Grammar == "NOT_PARSEABLE") %>%
  select(LipidName)

parsed_in <- parsed %>%
  filter(Original.Name %in% recognized)

lipidsigr_tbl <- LipidSigR::lipid_annotation(parsed_in)

annot_parseable <- lipidsigr_tbl %>%
  transmute(LipidName = feature, Class = class, Category = Category) %>%
  mutate(
    Class = str_squish(str_remove_all(as.character(Class), "\\[|\\]")),
    Class = if_else(is.na(Class) | Class == "", "Unassigned", Class),
    Category = str_squish(as.character(Category)),
    Category = str_remove(Category, "\\s*\\[[^\\]]+\\]\\s*$"),
    Category = str_squish(str_remove_all(Category, "\\[|\\]")),
    Category = case_when(
      Category %in% c("GP", "Glycerophospholipids") ~ "Glycerophospholipids",
      Category %in% c("GL", "Glycerolipids")        ~ "Glycerolipids",
      Category %in% c("SP", "Sphingolipids")        ~ "Sphingolipids",
      Category %in% c("ST", "Sterol lipids", "Sterol Lipids", "SL") ~ "Sterol lipids",
      Category %in% c("FA", "Fatty acyls")          ~ "Fatty acyls",
      is.na(Category) | Category == ""              ~ "Unassigned",
      TRUE ~ Category
    )
  )

manual_class_map <- tibble::tribble(
  ~LipidName,                       ~Class,         ~Category,
  "BA 24:1;O2;T/19:0;O",            "BA",           "Fatty acyls",
  "C18-Sphinganine",                "Sphinganine",  "Sphingolipids",
  "GPNAE 20:3",                     "GPNAE",        "Fatty acyls",
  "SE 24:1;O4/23:1",                "SE",           "Sterol lipids",
  "SE 24:2;O4/21:0;O",              "SE",           "Sterol lipids",
  "ST 27:0;O",                      "ST",           "Sphingolipids",
  "VAE 18:2",                       "VAE",          "Fatty acyls",
  "VAE 28:7",                       "VAE",          "Fatty acyls"
)

annot_manual <- not_parseable %>%
  mutate(
    LipidName = str_squish(as.character(LipidName)),
    Class = case_when(
      str_detect(LipidName, "^AHexCer\\s*\\(") ~ "AHexCer",
      TRUE ~ NA_character_
    )
  ) %>%
  left_join(manual_class_map, by = "LipidName", suffix = c("", ".map")) %>%
  mutate(
    Class = coalesce(Class, Class.map, "Unassigned"),
    Category = coalesce(
      Category,
      case_when(
        Class %in% c("AHexCer", "Sphinganine", "ST") ~ "Sphingolipids",
        Class %in% c("SE")                          ~ "Sterol lipids",
        Class %in% c("BA", "VAE", "GPNAE")          ~ "Fatty acyls",
        TRUE ~ "Unassigned"
      )
    )
  ) %>%
  select(LipidName, Class, Category)

lipid_annot <- bind_rows(annot_parseable, annot_manual) %>%
  distinct(LipidName, .keep_all = TRUE)

missing_ann <- setdiff(expr_final$LipidName, lipid_annot$LipidName)
if (length(missing_ann) > 0) {
  lipid_annot <- bind_rows(
    lipid_annot,
    tibble(LipidName = missing_ann, Class = "Unassigned", Category = "Unassigned")
  )
}

expr_final <- expr_final %>%
  left_join(lipid_annot, by = "LipidName") %>%
  select(any_of(info_cols), Class, Category, all_of(bio_cols))

# ------------------------------------------------------------------------------
# 9) Save output
# ------------------------------------------------------------------------------
write_delim(expr_final, out_file, delim = "\t")

cat("\n------------------------------------------------\n")
cat("Saved to:", out_file, "\n")
cat("Final matrix:", nrow(expr_final), "features x", length(bio_cols), "samples\n")
cat("------------------------------------------------\n")
