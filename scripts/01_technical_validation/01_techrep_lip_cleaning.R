# ==============================================================================
# Project: Plasma EV multi-omics – Lipidomics (TECH REPS)
# Script:  01_techrep_lip_cleaning.R
# Author:  Nicole Grinovero
# Date:    Feb 2026
#
# Aim:
#   Clean and curate MS-FLO lipidomics tables for tech-replicate method validation
#   (multi-omics plasma EV workflow), for both C18 POS and C18 NEG.
#
# Input:
#   data/techrep/lipidomics/C18_{pos,neg}/Area_plasma_EV_multiomics_techrep_C18_{pos,neg}.txt
#     (MS-DIAL export, pre MS-FLO)
#   data/techrep/lipidomics/C18_{pos,neg}/Area_plasma_EV_multiomics_techrep_C18_{pos,neg}_MSFLO_processed.txt
#     (MS-FLO processed)
#
# Output:
#   data/techrep/lipidomics/processed/techrep_lipids_cleaned_C18_pos.txt
#   data/techrep/lipidomics/processed/techrep_lipids_cleaned_C18_neg.txt
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
THR_QC_CV <- 0.30

# ---- POS paths ----
pos_dir <- here("data","techrep","lipidomics","C18_pos")
pos_pre_path  <- file.path(pos_dir, "Area_plasma_EV_multiomics_techrep_C18_pos.txt")            # raw pre MS-FLO
pos_expr_path <- file.path(pos_dir, "Area_plasma_EV_multiomics_techrep_C18_pos_MSFLO_processed.txt")  # MS-FLO processed

# ---- NEG paths ----
neg_dir <- here("data","techrep","lipidomics","C18_neg")
neg_pre_path  <- file.path(neg_dir, "Area_plasma_EV_multiomics_techrep_C18_neg.txt")            # raw pre MS-FLO
neg_expr_path <- file.path(neg_dir, "Area_plasma_EV_multiomics_techrep_C18_neg_MSFLO_processed.txt")  # MS-FLO processed

# outputs
out_dir <- here("data","techrep","lipidomics","processed")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_pos <- file.path(out_dir, "techrep_lipids_cleaned_C18_pos.txt")
out_neg <- file.path(out_dir, "techrep_lipids_cleaned_C18_neg.txt")

CONTAMINANTS_REGEX <- c(
  # MS-DIAL / generic
  "^RIKEN",
  "from\\s+Blank",
  "chemical",
  "reagent",
  
  # plasticizers / lab contaminants
  "phthalate",
  "siloxane",
  "polyethylene",
  "polypropylene",
  "polymer",
  "PEG",
  
  # fatty amide contaminants
  "docosenamide",
  "erucamide",
  "oleamide",
  
  # detergents / surfactants
  "dodecylbenzenesulfonic",
  "benzenesulfonic",
  "surfactant",
  "detergent",
  
  # antioxidants / additives from plastics
  "pentaerythritol",
  "hydrocinnamate",
  "hindered\\s*phenol",
  "butylated\\s*hydroxytoluene",
  "\\bBHT\\b",
  
  # xenobiotics / drugs
  "norethisterone"
)


info_cols <- c(
  "Alignment ID", "LipidName",
  "Average Mz", "Average Rt(min)", "Adduct type", "Ontology", "INCHIKEY", "SMILES"
)

# ------------------------------------------------------------------------------
# 0.1) Helper: recover/clean names from MS-FLO output (same as 03_lipidomics/01_lip_cleaning.R)
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
        is_merge & `MS/MS matched` == "True_False" ~ Name_Original_L,
        is_merge & `MS/MS matched` == "False_True" ~ Name_Original_R,
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
# 0.2) Helper: lipid key for safe deduplication (same as 03_lipidomics/01_lip_cleaning.R)
# ------------------------------------------------------------------------------
make_lipid_key <- function(x) {
  x <- str_squish(as.character(x))
  is_lyso <- str_detect(x, "^LP")
  x[is_lyso] <- str_replace(x[is_lyso], "/0:0$", "")
  str_squish(x)
}

# ==============================================================================
# POS
# ==============================================================================

# ------------------------------------------------------------------------------
# 1) Load POS (pre + processed)
# ------------------------------------------------------------------------------
pos_pre  <- read_delim(pos_pre_path,  delim = "\t", skip = 4, show_col_types = FALSE)
pos_expr <- read_delim(pos_expr_path, delim = "\t", show_col_types = FALSE)

cat("\n[POS] Processed features:", nrow(pos_expr), "\n")

# ------------------------------------------------------------------------------
# 2) Lipid name curation + contaminant removal (POS)
# ------------------------------------------------------------------------------
pos_names <- clean_lipid_names(pos_expr, pos_pre) %>%
  filter(!str_detect(
    LipidName,
    regex(paste(CONTAMINANTS_REGEX, collapse = "|"), ignore_case = TRUE)
  ))

cat("[POS] After name curation + contaminant removal:", nrow(pos_names), "\n")

# ------------------------------------------------------------------------------
# 3) Identify columns (bio/QC/BLK) and exclude ID_* (POS)
# ------------------------------------------------------------------------------
pos_sample_like <- setdiff(colnames(pos_names), info_cols)

pos_id_cols  <- pos_sample_like[str_detect(pos_sample_like, "_ID_")]
pos_qc_cols  <- pos_sample_like[str_detect(pos_sample_like, "_QC($|_|\\d)")]
pos_blk_cols <- pos_sample_like[str_detect(pos_sample_like, "_(BLK|BLANK)($|_)")]

pos_bio_cols <- pos_sample_like[
  str_detect(pos_sample_like, "^MagNet_") &
    !pos_sample_like %in% c(pos_id_cols, pos_qc_cols, pos_blk_cols)
]

cat("[POS] -> bio:", length(pos_bio_cols),
    "| QC:", length(pos_qc_cols),
    "| BLK:", length(pos_blk_cols),
    "| ID:", length(pos_id_cols), "\n")

# ------------------------------------------------------------------------------
# 4) Define missingness (non-positive -> NA) (POS)
# ------------------------------------------------------------------------------
pos_names <- pos_names %>%
  mutate(across(all_of(c(pos_bio_cols, pos_qc_cols, pos_blk_cols)), ~ suppressWarnings(as.numeric(.)))) %>%
  mutate(across(all_of(c(pos_bio_cols, pos_qc_cols, pos_blk_cols)), ~ if_else(. <= 0, NA_real_, .)))

# ------------------------------------------------------------------------------
# 5) QC filtering (CV over QC injections) (POS)
# ------------------------------------------------------------------------------
pos_qc <- pos_names %>%
  rowwise() %>%
  mutate(
    qc_mean = mean(c_across(all_of(pos_qc_cols)), na.rm = TRUE),
    qc_cv   = if_else(qc_mean > 0,
                      sd(c_across(all_of(pos_qc_cols)), na.rm = TRUE) / qc_mean,
                      NA_real_)
  ) %>%
  ungroup() %>%
  filter(!is.na(qc_cv) & qc_cv <= THR_QC_CV)

cat("[POS] After QC filtering (QC CV <=", THR_QC_CV, "):", nrow(pos_qc), "\n")

# ------------------------------------------------------------------------------
# 6) Duplicate resolution (one feature per LipidKey) (POS)
# ------------------------------------------------------------------------------
pos_qc2 <- pos_qc %>%
  mutate(LipidKey = make_lipid_key(LipidName))

pos_before <- nrow(pos_qc2)

pos_final <- pos_qc2 %>%
  group_by(LipidKey) %>%
  arrange(qc_cv, desc(qc_mean)) %>%   # no fill-rate criterion for technical replicates
  slice(1) %>%
  ungroup() %>%
  mutate(LipidName = LipidKey, Polarity = "C18_pos") %>%
  select(-LipidKey)

cat("[POS] Duplicate collapse:", pos_before, "->", nrow(pos_final), "features\n")

# ==============================================================================
# NEG
# ==============================================================================

# ------------------------------------------------------------------------------
# 1) Load NEG (pre + processed)
# ------------------------------------------------------------------------------
neg_pre  <- read_delim(neg_pre_path,  delim = "\t", skip = 4, show_col_types = FALSE)
neg_expr <- read_delim(
  neg_expr_path, delim = "\t", show_col_types = FALSE,
  col_types = cols(`MS/MS matched` = col_character())
)

cat("\n[NEG] Processed features:", nrow(neg_expr), "\n")

# ------------------------------------------------------------------------------
# 2) Lipid name curation + contaminant removal (NEG)
# ------------------------------------------------------------------------------
neg_names <- clean_lipid_names(neg_expr, neg_pre) %>%
  filter(!str_detect(
    LipidName,
    regex(paste(CONTAMINANTS_REGEX, collapse = "|"), ignore_case = TRUE)
  ))

cat("[NEG] After name curation + contaminant removal:", nrow(neg_names), "\n")

# ------------------------------------------------------------------------------
# 3) Identify columns (bio/QC/BLK) and exclude ID_* (NEG)
# ------------------------------------------------------------------------------
neg_sample_like <- setdiff(colnames(neg_names), info_cols)

neg_id_cols  <- neg_sample_like[str_detect(neg_sample_like, "_ID_")]
neg_qc_cols  <- neg_sample_like[str_detect(neg_sample_like, "_QC($|_|\\d)")]
neg_blk_cols <- neg_sample_like[str_detect(neg_sample_like, "_(BLK|BLANK)($|_)")]

neg_bio_cols <- neg_sample_like[
  str_detect(neg_sample_like, "^MagNet_") &
    !neg_sample_like %in% c(neg_id_cols, neg_qc_cols, neg_blk_cols)
]

cat("[NEG] -> bio:", length(neg_bio_cols),
    "| QC:", length(neg_qc_cols),
    "| BLK:", length(neg_blk_cols),
    "| ID:", length(neg_id_cols), "\n")

# ------------------------------------------------------------------------------
# 4) Define missingness (non-positive -> NA) (NEG)
# ------------------------------------------------------------------------------
neg_names <- neg_names %>%
  mutate(across(all_of(c(neg_bio_cols, neg_qc_cols, neg_blk_cols)), ~ suppressWarnings(as.numeric(.)))) %>%
  mutate(across(all_of(c(neg_bio_cols, neg_qc_cols, neg_blk_cols)), ~ if_else(. <= 0, NA_real_, .)))

# ------------------------------------------------------------------------------
# 5) QC filtering (CV over QC injections) (NEG)
# ------------------------------------------------------------------------------
neg_qc <- neg_names %>%
  rowwise() %>%
  mutate(
    qc_mean = mean(c_across(all_of(neg_qc_cols)), na.rm = TRUE),
    qc_cv   = if_else(qc_mean > 0,
                      sd(c_across(all_of(neg_qc_cols)), na.rm = TRUE) / qc_mean,
                      NA_real_)
  ) %>%
  ungroup() %>%
  filter(!is.na(qc_cv) & qc_cv <= THR_QC_CV)

cat("[NEG] After QC filtering (QC CV <=", THR_QC_CV, "):", nrow(neg_qc), "\n")

# ------------------------------------------------------------------------------
# 6) Duplicate resolution (one feature per LipidKey) (NEG)
# ------------------------------------------------------------------------------
neg_qc2 <- neg_qc %>%
  mutate(LipidKey = make_lipid_key(LipidName))

neg_before <- nrow(neg_qc2)

neg_final <- neg_qc2 %>%
  group_by(LipidKey) %>%
  arrange(qc_cv, desc(qc_mean)) %>%   # no fill-rate criterion for technical replicates
  slice(1) %>%
  ungroup() %>%
  mutate(LipidName = LipidKey, Polarity = "C18_neg") %>%
  select(-LipidKey)

cat("[NEG] Duplicate collapse:", neg_before, "->", nrow(neg_final), "features\n")

# ==============================================================================
# 7) Lipid class annotation (Goslin + LipidSigR) — POS
# ==============================================================================
parsed_pos <- rgoslin::parseLipidNames(lipidNames = pos_final$LipidName)

goslin_pos <- parsed_pos %>%
  transmute(LipidName = Original.Name, Grammar = Grammar)

recognized_pos <- goslin_pos %>% filter(Grammar != "NOT_PARSEABLE") %>% pull(LipidName)
not_parseable_pos <- goslin_pos %>% filter(Grammar == "NOT_PARSEABLE") %>% select(LipidName)

parsed_pos_in <- parsed_pos %>% filter(Original.Name %in% recognized_pos)
lipidsigr_pos <- LipidSigR::lipid_annotation(parsed_pos_in)

annot_pos_parseable <- lipidsigr_pos %>%
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

annot_pos_manual <- not_parseable_pos %>%
  mutate(
    LipidName = str_squish(as.character(LipidName)),
    Class = case_when(
      str_detect(LipidName, "^AHexCer\\s*\\(") ~ "AHexCer",
      str_detect(LipidName, "^SE\\s") | str_detect(LipidName, "^SE\\(") ~ "SE",
      TRUE ~ "Unassigned"
    ),
    Category = case_when(
      Class == "AHexCer" ~ "Sphingolipids",
      Class == "SE"      ~ "Sterol lipids",
      TRUE ~ "Unassigned"
    )
  ) %>%
  select(LipidName, Class, Category)

lip_annot_pos <- bind_rows(annot_pos_parseable, annot_pos_manual) %>%
  distinct(LipidName, .keep_all = TRUE)

missing_pos <- setdiff(pos_final$LipidName, lip_annot_pos$LipidName)
if (length(missing_pos) > 0) {
  lip_annot_pos <- bind_rows(
    lip_annot_pos,
    tibble(LipidName = missing_pos, Class = "Unassigned", Category = "Unassigned")
  )
}

pos_final <- pos_final %>%
  left_join(lip_annot_pos, by = "LipidName") %>%
  select(any_of(info_cols), Class, Category, Polarity, all_of(pos_bio_cols), all_of(pos_qc_cols), all_of(pos_blk_cols))

# ==============================================================================
# 8) Lipid class annotation (Goslin + LipidSigR) — NEG
# ==============================================================================
parsed_neg <- rgoslin::parseLipidNames(lipidNames = neg_final$LipidName)

goslin_neg <- parsed_neg %>%
  transmute(LipidName = Original.Name, Grammar = Grammar)

recognized_neg <- goslin_neg %>% filter(Grammar != "NOT_PARSEABLE") %>% pull(LipidName)
not_parseable_neg <- goslin_neg %>% filter(Grammar == "NOT_PARSEABLE") %>% select(LipidName)

parsed_neg_in <- parsed_neg %>% filter(Original.Name %in% recognized_neg)
lipidsigr_neg <- LipidSigR::lipid_annotation(parsed_neg_in)

annot_neg_parseable <- lipidsigr_neg %>%
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

annot_neg_manual <- not_parseable_neg %>%
  mutate(
    LipidName = str_squish(as.character(LipidName)),
    Class = case_when(
      str_detect(LipidName, "^AHexCer\\s*\\(") ~ "AHexCer",
      str_detect(LipidName, "^SE\\s") | str_detect(LipidName, "^SE\\(") ~ "SE",
      str_detect(LipidName, "^ST\\s") | str_detect(LipidName, "^ST\\(") ~ "ST",
      LipidName %in% c("Juniperic acid", "Ricinoleic acid") ~ "FA",
      str_detect(LipidName, "^NAGlySer\\b") ~ "NAGlySer",
      TRUE ~ "Unassigned"
    ),
    Category = case_when(
      Class == "AHexCer"  ~ "Sphingolipids",
      Class %in% c("SE","ST") ~ "Sterol lipids",
      Class %in% c("FA","NAGlySer") ~ "Fatty acyls",
      TRUE ~ "Unassigned"
    )
  ) %>%
  select(LipidName, Class, Category)

lip_annot_neg <- bind_rows(annot_neg_parseable, annot_neg_manual) %>%
  distinct(LipidName, .keep_all = TRUE)

missing_neg <- setdiff(neg_final$LipidName, lip_annot_neg$LipidName)
if (length(missing_neg) > 0) {
  lip_annot_neg <- bind_rows(
    lip_annot_neg,
    tibble(LipidName = missing_neg, Class = "Unassigned", Category = "Unassigned")
  )
}

neg_final <- neg_final %>%
  left_join(lip_annot_neg, by = "LipidName") %>%
  select(any_of(info_cols), Class, Category, Polarity, all_of(neg_bio_cols), all_of(neg_qc_cols), all_of(neg_blk_cols))

# ------------------------------------------------------------------------------
# 9) Save outputs
# ------------------------------------------------------------------------------
write_delim(pos_final, out_pos, delim = "\t")
write_delim(neg_final, out_neg, delim = "\t")

cat("\n------------------------------------------------\n")
cat("Saved POS to:", out_pos, "\n")
cat("Final POS:", nrow(pos_final), "features x", length(pos_bio_cols), "bio samples\n")
cat("Saved NEG to:", out_neg, "\n")
cat("Final NEG:", nrow(neg_final), "features x", length(neg_bio_cols), "bio samples\n")
cat("------------------------------------------------\n")
