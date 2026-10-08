# ==============================================================================
# Project: Preterm EV Multi-Omics – Proteomics
# Script:  02_prot_norm_impute.R
# Author:  Nicole Grinovero
# Date:    Jan 2026
#
# Aim:
#   From the cleaned protein matrix (Script 01), perform:
#     1) Biological filtering by timepoint (keep if >=70% observed in >=1 TP)
#     2) Sample-wise median realignment (log2 intensities, see 01_prot_cleaning.R)
#     3) Imputation in two steps:
#          - SLSA for MAR values only
#          - Downshifted Gaussian for MNAR values
#
# Inputs:
#   - data/proteomics/processed/proteins_cleaned.txt
#   - data/metadata/sample_metadata.txt
#
# Outputs:
#   - data/proteomics/processed/proteins_filtered.txt     (after filtering, before imputation; Fig. S7)
#   - data/proteomics/processed/mnar_mask.rds             (NA cells classified as MNAR)
#   - data/proteomics/processed/proteins_norm_imputed.txt
#
# Notes:
#   - MAR vs MNAR is inferred heuristically per timepoint block based on coverage.
#   - “MNAR” missingness is assumed left-censoring → imputed with downshifted N().
#   - QC plots are displayed in R, not saved.
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(imp4p)
})

set.seed(123)

# ------------------------------------------------------------------------------
# 0) Parameters
# ------------------------------------------------------------------------------
# (A) Biological filter: keep proteins with >=70% observed in at least one TP
THR_FILL <- 0.70

# (B) MNAR definition (timepoint-aware):
MNAR_N_CUTOFF <- 12
MNAR_THR_HIGH_N <- 0.30 # if n samples in TP >= 12
MNAR_THR_LOW_N <- 0.40 # if n samples in TP <  12

# (C) SLSA parameters:
SLSA_NKNN <- 30
SLSA_WEIGHT <- "o"
SLSA_INDCOMP <- 1

# (D) MNAR imputation (left-censored):
# Draw missing values from N(mean - IMP_SHIFT*sd, (IMP_SCALE*sd)^2)
IMP_SHIFT <- 1.8
IMP_SCALE <- 0.3

tp_levels <- c("T0", "T1", "T2", "T3", "T4")

# ------------------------------------------------------------------------------
# 1) Load data & metadata
# ------------------------------------------------------------------------------
path_clean <- here("data", "proteomics", "processed", "proteins_cleaned.txt")
path_meta <- here("data", "metadata", "sample_metadata.txt")

data <- read_delim(path_clean, delim = "\t", show_col_types = FALSE)
meta <- read_delim(path_meta, delim = "\t", show_col_types = FALSE)

cat("Input table:", nrow(data), "proteins x", ncol(data), "columns\n")

# ------------------------------------------------------------------------------
# 2) Identify intensity columns (biological sample_id) and align metadata
# ------------------------------------------------------------------------------
# After Script 01, sample columns are renamed to biological sample_id.
meta_use <- meta %>%
  mutate(
    sample_id = as.character(sample_id),
    timepoint = as.character(timepoint)
  ) %>%
  filter(sample_id %in% names(data))

if (nrow(meta_use) == 0) {
  stop("No sample_id from metadata found in proteins_cleaned.txt. Check Script 01 output / metadata.")
}

# Prefer explicit T0–T4 ordering if applicable
if (all(unique(meta_use$timepoint) %in% tp_levels)) {
  meta_use <- meta_use %>% mutate(timepoint = factor(timepoint, levels = tp_levels))
} else {
  warning("Timepoints in metadata are not a subset of {T0..T4}. Using alphabetical factor order.")
  meta_use <- meta_use %>% mutate(timepoint = factor(timepoint))
}

# Define intensity columns in the SAME order as metadata (stable & reproducible)
meta_use <- meta_use %>% arrange(timepoint, sample_id)
intensity_cols <- meta_use$sample_id
timepoints <- levels(meta_use$timepoint)

cat("Samples detected:", length(intensity_cols), "\n")
cat("Timepoints detected:", paste(timepoints, collapse = ", "), "\n")

# ------------------------------------------------------------------------------
# 3) Biological filtering (timepoint-aware fill rate)
# ------------------------------------------------------------------------------
# Rationale:
# Keep proteins that are reasonably quantified in at least one timepoint,
# removing extremely sparse features that can destabilize modeling/clustering.
keep_row <- rep(FALSE, nrow(data))

for (tp in timepoints) {
  tp_samples <- meta %>%
    filter(timepoint == tp) %>%
    pull(sample_id)
  tp_cols <- intersect(intensity_cols, tp_samples)
  if (length(tp_cols) == 0) next

  sub_mat <- as.matrix(data[, tp_cols, drop = FALSE])
  prop_non_na <- rowMeans(!is.na(sub_mat))
  keep_row <- keep_row | (prop_non_na >= THR_FILL)
}

filtered_data <- data[keep_row, , drop = FALSE]
cat(
  "Proteins kept (fill-rate >= ", THR_FILL, " in >=1 TP): ",
  nrow(filtered_data), "\n", sep = ""
)

# Post-filter, pre-imputation table (input to Fig. S7)
write_delim(filtered_data,
            here("data", "proteomics", "processed", "proteins_filtered.txt"),
            delim = "\t")

# ------------------------------------------------------------------------------
# 4) Expression matrix (log2 intensities)
# ------------------------------------------------------------------------------
expr_matrix <- filtered_data %>%
  select(all_of(intensity_cols)) %>%
  as.matrix()

# Use UniProt as rownames if available (helps downstream tracking)
if ("UniProtID" %in% names(filtered_data)) {
  rownames(expr_matrix) <- filtered_data$UniProtID
}

na0 <- sum(is.na(expr_matrix))
cat(
  "Missing values after filtering:", na0,
  sprintf("(%.2f%%)\n", 100 * na0 / length(expr_matrix))
)

# QC: % missing values per timepoint
missing_df <- as.data.frame(expr_matrix) %>%
  pivot_longer(cols = everything(), names_to = "sample_id", values_to = "Intensity") %>%
  left_join(meta_use %>% select(sample_id, timepoint), by = "sample_id") %>%
  group_by(timepoint) %>%
  summarise(pct_missing = mean(is.na(Intensity)) * 100, .groups = "drop")

p_missing <- ggplot(missing_df, aes(x = timepoint, y = pct_missing, fill = timepoint)) +
  geom_col(width = 0.7, alpha = 0.9) +
  scale_fill_viridis_d(option = "D", end = 0.85, guide = "none") +
  labs(title = "Missingness by timepoint", x = "Timepoint", y = "% missing values") +
  theme_bw(base_size = 12) +
  ylim(0, 100)

print(p_missing)

# ------------------------------------------------------------------------------
# 5) Normalization: sample-wise median realignment
# ------------------------------------------------------------------------------
# Rationale:
# Align sample medians to reduce technical shifts while preserving
# relative differences between proteins.
sample_medians <- apply(expr_matrix, 2, median, na.rm = TRUE)
global_median  <- median(sample_medians, na.rm = TRUE)

expr_median_norm <- sweep(expr_matrix, 2, sample_medians, FUN = "-")
expr_median_norm <- expr_median_norm + global_median

# QC: pre- vs post-normalization boxplots
plot_data <- bind_rows(
  as_tibble(expr_matrix, rownames = "feature") %>%
    pivot_longer(-feature, names_to = "sample_id", values_to = "Intensity") %>%
    mutate(Status = "Raw"),
  as_tibble(expr_median_norm, rownames = "feature") %>%
    pivot_longer(-feature, names_to = "sample_id", values_to = "Intensity") %>%
    mutate(Status = "Median-aligned")
) %>%
  left_join(meta_use %>% select(sample_id, timepoint), by = "sample_id")

p_box <- ggplot(plot_data, aes(x = sample_id, y = Intensity, fill = timepoint)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.25) +
  facet_wrap(~ Status, ncol = 1, scales = "free_y") +
  scale_fill_viridis_d(option = "D", end = 0.9) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right"
  ) +
  labs(title = "Proteomics: raw vs normalized", x = NULL, y = "log2 intensity")

print(p_box)

# ------------------------------------------------------------------------------
# 6) MNAR mask (timepoint-aware)
# ------------------------------------------------------------------------------
# We classify a missing cell as MNAR if:
#   - the value is NA
#   - within that timepoint block, the protein has very low coverage 
mask_MNAR_tp <- matrix(FALSE,
                       nrow = nrow(expr_median_norm),
                       ncol = ncol(expr_median_norm),
                       dimnames = list(rownames(expr_median_norm), colnames(expr_median_norm))
)

for (tp in timepoints) {
  tp_cols <- meta_use %>% filter(timepoint == tp) %>% pull(sample_id)
  tp_cols <- intersect(tp_cols, colnames(expr_median_norm))
  if (length(tp_cols) == 0) next
  
  block <- expr_median_norm[, tp_cols, drop = FALSE]
  prop_non_na <- rowMeans(!is.na(block))
  
  thr <- if (length(tp_cols) >= MNAR_N_CUTOFF) MNAR_THR_HIGH_N else MNAR_THR_LOW_N
  mnar_rows <- prop_non_na < thr
  if (!any(mnar_rows)) next
  
  mask_MNAR_tp[mnar_rows, tp_cols] <- is.na(block)[mnar_rows, , drop = FALSE]
}

cat("NA cells classified as MNAR:", sum(mask_MNAR_tp), "\n")
saveRDS(mask_MNAR_tp, file = here("data", "proteomics", "processed", "mnar_mask.rds"))


# ------------------------------------------------------------------------------
# 7) Imputation part 1: SLSA for MAR only
# ------------------------------------------------------------------------------
# Strategy:
# - Run SLSA on the normalized matrix.
# - Copy imputed values ONLY into NA cells that are NOT classified as MNAR.
conditions <- meta_use %>%
  arrange(match(sample_id, colnames(expr_median_norm))) %>%
  pull(timepoint)

expr_slsa_full <- impute.slsa(
  tab          = expr_median_norm,
  conditions   = conditions,
  nknn         = SLSA_NKNN,
  selec        = "all",
  weight       = SLSA_WEIGHT,
  ind.comp     = SLSA_INDCOMP,
  progress.bar = TRUE
)

expr_after_slsa <- expr_median_norm
to_impute_MAR <- is.na(expr_after_slsa) & !mask_MNAR_tp
cat("Cells imputed with SLSA (MAR):", sum(to_impute_MAR), "\n")

expr_after_slsa[to_impute_MAR] <- expr_slsa_full[to_impute_MAR]

# ------------------------------------------------------------------------------
# 8) Imputation part 2: downshifted Gaussian for MNAR
# ------------------------------------------------------------------------------
expr_final <- expr_after_slsa
to_impute_MNAR <- is.na(expr_final) & mask_MNAR_tp
cat("Cells imputed with downshifted Gaussian (MNAR):", sum(to_impute_MNAR), "\n")

for (j in seq_len(ncol(expr_final))) {
  x <- expr_final[, j]
  mnar_j <- is.na(x) & mask_MNAR_tp[, j]
  if (!any(mnar_j)) next
  
  valid_vals <- x[!is.na(x)]
  if (length(valid_vals) < 5) {
    # extremely rare after SLSA, but keeps behavior safe
    next
  }
  
  mu <- mean(valid_vals)
  sigma <- sd(valid_vals)
  if (is.na(sigma) || sigma == 0) sigma <- 1e-6
  
  mu_imp <- mu - (IMP_SHIFT * sigma)
  sd_imp <- sigma * IMP_SCALE
  
  x[mnar_j] <- rnorm(sum(mnar_j), mean = mu_imp, sd = sd_imp)
  expr_final[, j] <- x
}

stopifnot(sum(is.na(expr_final)) == 0)

# ------------------------------------------------------------------------------
# 9) QC: PCA (post-normalization + imputation)
# ------------------------------------------------------------------------------
pca_res <- prcomp(t(expr_final), scale. = FALSE)

pca_plot_df <- as_tibble(pca_res$x, rownames = "sample_id") %>%
  left_join(meta_use, by = "sample_id")

var_expl <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)

p_pca <- ggplot(pca_plot_df, aes(x = PC1, y = PC2, color = timepoint)) +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_viridis_d(option = "D", end = 0.9) +
  labs(
    title = "Proteomics – PCA (post norm + imputation)",
    subtitle = "Median realignment + SLSA (MAR) + Downshifted Gaussian (MNAR)",
    x = paste0("PC1 (", var_expl[1], "%)"),
    y = paste0("PC2 (", var_expl[2], "%)")
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank())

print(p_pca)

# ------------------------------------------------------------------------------
# 10) Reassemble final dataframe and save
# ------------------------------------------------------------------------------
info_cols <- setdiff(names(filtered_data), intensity_cols)

final_prot_df <- bind_cols(
  filtered_data %>% select(all_of(info_cols)),
  as_tibble(expr_final)
)

out_dir <- here("data", "proteomics", "processed")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_file <- file.path(out_dir, "proteins_norm_imputed.txt")
write_delim(final_prot_df, out_file, delim = "\t")

cat("Saved to:", out_file, "\n")
