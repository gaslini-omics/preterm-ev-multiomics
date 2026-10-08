# ==============================================================================
# Project: Preterm EV Multi-Omics – Proteomics
# Script:  04_prot_stats_time_LMM.R
# Author:  Nicole Grinovero
# Date:    Jan 2026
#
# Aim:
#   Identify EV-associated proteins that change from birth (T0)
#   to term-equivalent age (T4) using a linear mixed model (LMM).
#
# Model (per protein):
#   intensity ~ timepoint + ga_birth + (1 | patient_id)
#
# Tests:
#   - Global time effect: Type-III ANOVA on "timepoint"
#   - Primary endpoint: coefficient for timepointT4 (T4 vs T0 reference)
#
# Output:
#   outputs/tables/preterm/proteomics/proteins_LMM_time_results.txt
#
# Notes:
#   - timepoint is treated as a categorical factor with T0 as reference level.
#   - patient_id random intercept accounts for repeated measures.
#   - ga_birth included as covariate (numeric).
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(lmerTest)
})

# ------------------------------------------------------------------------------
# 0) Configuration
# ------------------------------------------------------------------------------
FDR_CUTOFF <- 0.05
tp_levels <- c("T0", "T1", "T2", "T3", "T4")

expr_path <- here("data", "proteomics", "processed", "proteins_norm_imputed.txt")
meta_path <- here("data", "metadata", "sample_metadata.txt")

out_dir <- here("outputs", "tables", "preterm", "proteomics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, "proteins_LMM_time_results.txt")

ID_COL <- "UniProtID"

# ------------------------------------------------------------------------------
# 1) Data preparation
# ------------------------------------------------------------------------------
expr <- read_delim(expr_path, delim = "\t", show_col_types = FALSE)
meta <- read_delim(meta_path, delim = "\t", show_col_types = FALSE)

# Required columns
req_meta <- c("sample_id", "timepoint", "patient_id", "ga_birth")
req_expr <- c(ID_COL, "GeneName")

stopifnot(all(req_meta %in% names(meta)))
stopifnot(all(req_expr %in% names(expr)))

# Detect sample columns in expression matrix (preserve metadata order!)
bio_cols <- meta$sample_id[meta$sample_id %in% names(expr)]
if (length(bio_cols) == 0) {
  stop("No metadata sample_id columns found in expression table. Check inputs.")
}

info_cols <- setdiff(names(expr), bio_cols)

# Extract annotation
annot <- expr %>%
  select(all_of(info_cols)) %>%
  distinct(.data[[ID_COL]], .keep_all = TRUE)

# Reshape to long format for modeling
expr_long <- expr %>%
  pivot_longer(all_of(bio_cols), names_to = "sample_id", values_to = "intensity") %>%
  left_join(meta, by = "sample_id") %>%
  mutate(
    timepoint  = factor(timepoint, levels = c("T0", "T1", "T2", "T3", "T4")),
    patient_id = factor(patient_id),
    ga_birth   = as.numeric(ga_birth)
  )

cat("Proteins (unique IDs):", n_distinct(expr_long[[ID_COL]]), "\n")
cat("Observations (rows):", nrow(expr_long), "\n")

# ------------------------------------------------------------------------------
# 2) Statistical modeling (LMM per protein)
# ------------------------------------------------------------------------------
# We fit the LMM for each protein and extract:
#  - p_global: ANOVA p-value for timepoint
#  - beta_T4_vs_T0 and p_T4_vs_T0: coefficient for timepointT4 (requires T0 reference)
fit_lmm <- function(sub_df) {
  ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))

  tryCatch(
    {
      m <- lmer(intensity ~ timepoint + ga_birth + (1 | patient_id),
        data = sub_df, REML = FALSE, control = ctrl
      )

      # Extract stats
      aov_res <- anova(m, type = 3)
      coefs <- summary(m)$coefficients

      tibble(
        p_global      = aov_res["timepoint", "Pr(>F)"],
        beta_T4_vs_T0 = coefs["timepointT4", "Estimate"],
        p_T4_vs_T0    = coefs["timepointT4", "Pr(>|t|)"]
      )
    },
    error = function(e) tibble(p_global = NA, beta_T4_vs_T0 = NA, p_T4_vs_T0 = NA)
  )
}

# Run loop using purrr::map (split by UniProtID)
results_raw <- split(expr_long, expr_long[[ID_COL]]) %>%
  map_dfr(fit_lmm, .id = ID_COL)

# ------------------------------------------------------------------------------
# 3) Multiple testing correction + results export
# ------------------------------------------------------------------------------
results_final <- results_raw %>%
  left_join(annot, by = ID_COL) %>%
  mutate(
    FDR_global = p.adjust(p_global, method = "BH"),
    FDR_T4_T0  = p.adjust(p_T4_vs_T0, method = "BH")
  ) %>%
  relocate(
    UniProtID,
    GeneName,
    beta_T4_vs_T0,
    p_T4_vs_T0,
    FDR_T4_T0,
    p_global,
    FDR_global,
    .before = 1
  ) %>%
  arrange(FDR_global)

# Save results
write_delim(results_final, out_file, delim = "\t")

# Summaries
n_total <- nrow(results_final)
n_sig_global <- sum(results_final$FDR_global < FDR_CUTOFF, na.rm = TRUE)
n_sig_t4 <- sum(results_final$FDR_T4_T0 < FDR_CUTOFF, na.rm = TRUE)

cat("\n------------------------------------------------\n")
cat(">>> PROTEOMICS ANALYSIS COMPLETE. Summary:\n")
cat("    Total proteins analyzed :", n_total, "\n")
cat("    Global Time Effect      :", n_sig_global, "(FDR < 0.05)\n")
cat("    Changed at T4 vs T0     :", n_sig_t4, "(FDR < 0.05)\n")
cat("------------------------------------------------\n")
cat("Results saved to:", out_file, "\n")
