# ==============================================================================
# Project: Preterm EV Multi-Omics – Lipidomics
# Script:  03_lip_stats_time_LMM.R
# Author:  Nicole Grinovero
# Date:    Jan 2026
#
# Aim:
#   Identify EV-associated lipids that change from birth (T0) to term-equivalent
#   age (T4) using a linear mixed model (LMM).
#
# Model (per lipid):
#   intensity ~ timepoint + ga_birth + (1 | patient_id)
#
# Tests:
#   - Global time effect: Type-III ANOVA on "timepoint"
#   - Primary endpoint: coefficient for timepointT4 (T4 vs T0 reference)
#
# Output:
#   outputs/tables/preterm/lipidomics/lipids_LMM_time_results.txt
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(lmerTest)
})

# ------------------------------------------------------------------------------
# 0) Configuration
# ------------------------------------------------------------------------------
FDR_CUTOFF <- 0.05
tp_levels  <- c("T0", "T1", "T2", "T3", "T4")

expr_path <- here("data", "lipidomics", "processed", "lipids_norm_imputed.txt")
meta_path <- here("data", "metadata", "sample_metadata.txt")

out_dir  <- here("outputs", "tables", "preterm", "lipidomics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, "lipids_LMM_time_results.txt")

ID_COL <- "LipidName"

# ------------------------------------------------------------------------------
# 1) Data preparation
# ------------------------------------------------------------------------------
# Deuterated internal standards ("(dX)" in the name) are excluded
expr <- read_delim(expr_path, delim = "\t", show_col_types = FALSE) %>% 
  filter(!str_detect(LipidName, "\\(d[0-9]+\\)"))

meta <- read_delim(meta_path, delim = "\t", show_col_types = FALSE)

req_meta <- c("sample_id", "timepoint", "patient_id", "ga_birth")
stopifnot(all(req_meta %in% names(meta)))
stopifnot(ID_COL %in% names(expr))

# Detect sample columns in metadata order (reproducible)
bio_cols <- meta$sample_id[meta$sample_id %in% names(expr)]
if (length(bio_cols) == 0) {
  stop("No metadata sample_id columns found in lipid expression table. Check inputs.")
}

info_cols <- setdiff(names(expr), bio_cols)

# LipidName should be unique after cleaning; stop early if not
dup_ids <- expr[[ID_COL]][duplicated(expr[[ID_COL]])]
if (length(dup_ids) > 0) {
  stop("Duplicate LipidName values found (n=", length(unique(dup_ids)),
       "). Resolve before LMM (dedup step in 01_lip_cleaning).")
}

# Annotation (one row per LipidName)
annot <- expr %>%
  select(all_of(info_cols)) %>%
  distinct(.data[[ID_COL]], .keep_all = TRUE)

# Long format for LMM
expr_long <- expr %>%
  pivot_longer(all_of(bio_cols), names_to = "sample_id", values_to = "intensity") %>%
  left_join(meta, by = "sample_id") %>%
  mutate(
    timepoint  = factor(as.character(timepoint), levels = tp_levels),
    patient_id = factor(patient_id),
    ga_birth   = as.numeric(ga_birth)
  )

cat("Lipids (unique):", n_distinct(expr_long[[ID_COL]]), "\n")
cat("Observations:", nrow(expr_long), "\n")

# ------------------------------------------------------------------------------
# 2) Fit LMM per lipid
# ------------------------------------------------------------------------------
fit_lmm <- function(sub_df) {
  ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
  
  tryCatch(
    {
      m <- lmer(
        intensity ~ timepoint + ga_birth + (1 | patient_id),
        data = sub_df, REML = FALSE, control = ctrl
      )
      
      aov_res <- anova(m, type = 3)
      coefs   <- summary(m)$coefficients
      
      tibble(
        p_global      = aov_res["timepoint", "Pr(>F)"],
        beta_T4_vs_T0 = coefs["timepointT4", "Estimate"],
        p_T4_vs_T0    = coefs["timepointT4", "Pr(>|t|)"]
      )
    },
    error = function(e) tibble(p_global = NA_real_, beta_T4_vs_T0 = NA_real_, p_T4_vs_T0 = NA_real_)
  )
}

results_raw <- split(expr_long, expr_long[[ID_COL]]) %>%
  map_dfr(fit_lmm, .id = ID_COL)

# ------------------------------------------------------------------------------
# 3) Multiple testing correction + results export
# ------------------------------------------------------------------------------
results_final <- results_raw %>%
  mutate(
    FDR_global = p.adjust(p_global, method = "BH"),
    FDR_T4_T0  = p.adjust(p_T4_vs_T0, method = "BH")
  ) %>%
  left_join(annot, by = ID_COL) %>%
  relocate(
    # IDs
    `Alignment ID`,
    LipidName,
    # primary endpoint
    beta_T4_vs_T0,
    p_T4_vs_T0,
    FDR_T4_T0,
    # global test
    p_global,
    FDR_global,
    .before = 1
  ) %>%
  arrange(FDR_global)

write_delim(results_final, out_file, delim = "\t")

# ------------------------------------------------------------------------------
# 4) Summary
# ------------------------------------------------------------------------------
n_total      <- nrow(results_final)
n_sig_global <- sum(results_final$FDR_global < FDR_CUTOFF, na.rm = TRUE)
n_sig_t4     <- sum(results_final$FDR_T4_T0 < FDR_CUTOFF, na.rm = TRUE)

cat("\n------------------------------------------------\n")
cat(">>> LIPIDOMICS ANALYSIS COMPLETE. Summary:\n")
cat("    Total lipids analyzed :", n_total, "\n")
cat("    Global Time Effect    :", n_sig_global, "(FDR < 0.05)\n")
cat("    Changed at T4 vs T0   :", n_sig_t4, "(FDR < 0.05)\n")
cat("------------------------------------------------\n")
cat("Results saved to:", out_file, "\n")
