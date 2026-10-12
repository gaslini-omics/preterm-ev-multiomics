# ==============================================================================
# Project: Preterm EV Multi-Omics – Inputs for per-timepoint comparisons
# Script:  04_export_statistics_inputs.R
# Author:  Nicole Grinovero
#
# Aim:
#   Export the input tables used by scripts/06_statistics/ (per-timepoint
#   Mann-Whitney comparisons, Supplementary Table S8): the leading-edge features
#   of the significant module–outcome associations (FDR_global < 0.05), their
#   normalized and imputed abundances, and the binary clinical outcomes.
#
# Inputs:
#   - outputs/tables/preterm/integration/clinical_camera/modules_camera_all_outcomes_multiomics.tsv
#   - outputs/tables/preterm/integration/clinical_camera/drivers_top2_<layer>_<outcome>.tsv
#   - data/proteomics/processed/proteins_norm_imputed.txt
#   - data/lipidomics/processed/lipids_norm_imputed.txt
#   - data/metadata/sample_metadata.txt
#
# Outputs (outputs/tables/preterm/statistics/input/):
#   - drivers_summary.tsv          Layer, ModuleLabel, Outcome, FeatureName
#   - drivers_expr_proteins.tsv    UniProtID, FeatureName (= gene name), one column per sample
#   - drivers_expr_lipids.tsv      FeatureName (= lipid name), one column per sample
#   - outcomes_by_sample.tsv       sample_id, timepoint, binary clinical outcomes
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(purrr)
})

# ------------------------------------------------------------------------------
# 0) Configuration
# ------------------------------------------------------------------------------
OUTCOMES <- c(
  "any_brain_lesion",
  "hemorrhagic_lesion",
  "sepsis",
  "mechanical_ventilation"
)
LAYERS     <- c("Protein", "Lipid")
FDR_CUTOFF <- 0.05   # same threshold used for Fig. 7a

cam_dir   <- here("outputs", "tables", "preterm", "integration", "clinical_camera")
in_master <- file.path(cam_dir, "modules_camera_all_outcomes_multiomics.tsv")

prot_expr_path <- here("data", "proteomics", "processed", "proteins_norm_imputed.txt")
lip_expr_path  <- here("data", "lipidomics",  "processed", "lipids_norm_imputed.txt")
meta_path      <- here("data", "metadata", "sample_metadata.txt")

out_dir <- here("outputs", "tables", "preterm", "statistics", "input")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1) Significant module–outcome associations
# ------------------------------------------------------------------------------
sel <- read_tsv(in_master, show_col_types = FALSE) %>%
  mutate(FDR_global = as.numeric(FDR_global)) %>%
  filter(FDR_global < FDR_CUTOFF) %>%
  distinct(Layer, Outcome, ModuleLabel)

stopifnot(nrow(sel) > 0)
cat("Significant module–outcome associations (FDR_global <", FDR_CUTOFF, "):",
    nrow(sel), "\n")

# ------------------------------------------------------------------------------
# 2) Leading-edge features of those modules (drivers_top2 files from 03)
# ------------------------------------------------------------------------------
read_drivers <- function(layer, outcome) {
  f <- file.path(cam_dir, paste0("drivers_top2_", layer, "_", outcome, ".tsv"))
  stopifnot(file.exists(f))
  d <- read_tsv(f, show_col_types = FALSE)
  # proteins are identified by UniProt accession and reported by gene name
  d$FeatureName <- if (layer == "Protein") d$GeneName else d$Feature
  d %>% select(Layer, ModuleLabel, Outcome, Feature, FeatureName)
}

drivers <- sel %>%
  distinct(Layer, Outcome) %>%
  arrange(factor(Layer, levels = LAYERS), factor(Outcome, levels = OUTCOMES)) %>%
  pmap_dfr(\(Layer, Outcome) read_drivers(Layer, Outcome)) %>%
  semi_join(sel, by = c("Layer", "Outcome", "ModuleLabel"))

# every significant module must have its leading-edge features
stopifnot(nrow(anti_join(sel, drivers, by = c("Layer", "Outcome", "ModuleLabel"))) == 0)
stopifnot(!anyNA(drivers$FeatureName))

# ------------------------------------------------------------------------------
# 3) Abundances of the leading-edge features (normalized, imputed, log2)
# ------------------------------------------------------------------------------
meta <- read_tsv(meta_path, show_col_types = FALSE)

prot <- read_tsv(prot_expr_path, show_col_types = FALSE)
prot_samples <- names(prot)[names(prot) %in% meta$sample_id]
prot_out <- prot %>%
  filter(UniProtID %in% drivers$Feature[drivers$Layer == "Protein"]) %>%
  select(UniProtID, FeatureName = GeneName, all_of(prot_samples))

lip <- read_tsv(lip_expr_path, show_col_types = FALSE)
lip_samples <- names(lip)[names(lip) %in% meta$sample_id]
lip_out <- lip %>%
  filter(LipidName %in% drivers$Feature[drivers$Layer == "Lipid"]) %>%
  select(FeatureName = LipidName, all_of(lip_samples))

# one row per feature, and every driver found in the abundance tables
stopifnot(!anyDuplicated(prot_out$FeatureName), !anyDuplicated(lip_out$FeatureName))
stopifnot(setequal(prot_out$FeatureName, drivers$FeatureName[drivers$Layer == "Protein"]))
stopifnot(setequal(lip_out$FeatureName,  drivers$FeatureName[drivers$Layer == "Lipid"]))

# ------------------------------------------------------------------------------
# 4) Clinical outcomes by sample
# ------------------------------------------------------------------------------
outcomes_out <- meta %>%
  select(sample_id, timepoint, all_of(OUTCOMES))

# ------------------------------------------------------------------------------
# 5) Write
# ------------------------------------------------------------------------------
write_tsv(drivers %>% select(Layer, ModuleLabel, Outcome, FeatureName),
          file.path(out_dir, "drivers_summary.tsv"))
write_tsv(prot_out,     file.path(out_dir, "drivers_expr_proteins.tsv"))
write_tsv(lip_out,      file.path(out_dir, "drivers_expr_lipids.tsv"))
write_tsv(outcomes_out, file.path(out_dir, "outcomes_by_sample.tsv"))

cat("Leading-edge features:", nrow(drivers), "rows (",
    nrow(prot_out), "proteins,", nrow(lip_out), "lipids )\n")
cat("Saved in:", out_dir, "\n")
