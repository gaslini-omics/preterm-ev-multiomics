# ==============================================================================
# Project: Preterm EV Multi-Omics – Lipidomics
# Script:  02_lip_norm_impute.R
# Author:  Nicole Grinovero
# Date:    Jan 2026
#
# Aim:
#   Starting from the cleaned lipidomics matrix (Script 01), perform:
#     1) log2 transformation (areas -> log2 intensities)
#     2) sample-wise median normalization (correct global intensity shifts)
#     3) missing value imputation assuming left-censoring (MNAR; LC–MS typical)
#     4) QC visualizations (boxplots + PCA; shown in R, not saved)
#
# Notes:
#   - Missing values represent non-detections (below LOD/LOQ): MNAR assumption.
#   - Imputation is performed per sample by drawing values from a downshifted
#     Gaussian distribution: N(mu - shift*sigma, (scale*sigma)^2).
#
# Output:
#   data/lipidomics/processed/lipids_norm_imputed.txt
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

set.seed(123)

# ------------------------------------------------------------------------------
# 0) Configuration
# ------------------------------------------------------------------------------
tp_levels <- c("T0", "T1", "T2", "T3", "T4")

# MNAR / left-censored imputation parameters 
IMP_SHIFT <- 1.8   # mean downshift in SD units
IMP_SCALE <- 0.3   # SD scaling factor (narrower distribution)

proc_dir <- here("data", "lipidomics", "processed")
meta_path <- here("data", "metadata", "sample_metadata.txt")

in_file  <- file.path(proc_dir, "lipids_cleaned.txt")
out_file <- file.path(proc_dir, "lipids_norm_imputed.txt")

# ------------------------------------------------------------------------------
# 1) Load cleaned data + metadata
# ------------------------------------------------------------------------------
expr_clean <- read_delim(in_file, delim = "\t", show_col_types = FALSE)
meta <- read_delim(meta_path, delim = "\t", show_col_types = FALSE) %>%
  mutate(timepoint = factor(timepoint, levels = tp_levels))

stopifnot("sample_id" %in% names(meta))
stopifnot("timepoint" %in% names(meta))

# Use metadata order 
bio_cols <- meta$sample_id[meta$sample_id %in% colnames(expr_clean)]
if (length(bio_cols) == 0) {
  stop("No sample columns matched between lipids_cleaned and meta$sample_id.")
}

info_cols <- setdiff(colnames(expr_clean), bio_cols)

meta_lip <- meta %>%
  filter(sample_id %in% bio_cols) %>%
  arrange(match(sample_id, bio_cols))

# Ensure sample column order matches metadata
stopifnot(identical(meta_lip$sample_id, bio_cols))

# Extract numeric matrix (features x samples)
mat_raw <- expr_clean %>%
  select(all_of(bio_cols)) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

cat("Input matrix:", nrow(mat_raw), "features x", ncol(mat_raw), "samples\n")
na0 <- sum(is.na(mat_raw))
cat("Missing values:", na0, sprintf("(%.2f%%)\n", 100 * na0 / length(mat_raw)))

# ------------------------------------------------------------------------------
# 2) Log2 transformation
# ------------------------------------------------------------------------------
# Non-positive values are treated as missing
mat_raw[mat_raw <= 0] <- NA_real_
mat_log2 <- log2(mat_raw)

# ------------------------------------------------------------------------------
# 3) Sample-wise median normalization (median alignment)
# ------------------------------------------------------------------------------
# Rationale:
# Correct global shifts between samples (injection amount, batch intensity drift)
# while preserving relative differences between lipids within each sample.
sample_medians <- apply(mat_log2, 2, median, na.rm = TRUE)
global_median  <- median(sample_medians, na.rm = TRUE)

mat_norm <- sweep(mat_log2, 2, sample_medians, FUN = "-")
mat_norm <- mat_norm + global_median

# ------------------------------------------------------------------------------
# 4) QC plot: boxplots before vs after normalization (shown in R, not saved)
# ------------------------------------------------------------------------------
plot_data <- bind_rows(
  as_tibble(mat_log2, rownames = "feature") %>%
    pivot_longer(-feature, names_to = "sample_id", values_to = "Intensity") %>%
    mutate(Status = "1. Raw (log2)"),
  as_tibble(mat_norm, rownames = "feature") %>%
    pivot_longer(-feature, names_to = "sample_id", values_to = "Intensity") %>%
    mutate(Status = "2. Median-aligned")
) %>%
  left_join(meta_lip %>% select(sample_id, timepoint), by = "sample_id")

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
  labs(title = "Lipidomics: raw vs normalized", x = NULL, y = "log2 intensity")

print(p_box)

# ------------------------------------------------------------------------------
# 5) Missing value imputation (left-censored Gaussian, per sample)
# ------------------------------------------------------------------------------
# Rationale:
# In untargeted LC–MS, missingness is often due to low abundance (below detection).
# We therefore impute low values rather than mean-like values (MNAR assumption).
impute_left_censored <- function(x, shift = 1.8, scale = 0.3) {
  miss <- is.na(x)
  n_missing <- sum(miss)
  if (n_missing == 0) return(x)
  
  valid_vals <- x[!miss]
  
  # If too few observed values, fall back to robust defaults
  if (length(valid_vals) < 3) {
    mu <- mean(valid_vals, na.rm = TRUE)
    if (is.na(mu)) mu <- 0
    sigma <- sd(valid_vals, na.rm = TRUE)
    if (is.na(sigma) || sigma == 0) sigma <- 1e-6
  } else {
    mu <- mean(valid_vals)
    sigma <- sd(valid_vals)
    if (is.na(sigma) || sigma == 0) sigma <- 1e-6
  }
  
  mu_imp <- mu - (shift * sigma)
  sd_imp <- sigma * scale
  if (is.na(sd_imp) || sd_imp <= 0) sd_imp <- 1e-6
  
  x[miss] <- rnorm(n_missing, mean = mu_imp, sd = sd_imp)
  x
}

mat_imputed <- apply(
  mat_norm, 2,
  impute_left_censored,
  shift = IMP_SHIFT, scale = IMP_SCALE
)

# Re-attach dimnames 
dimnames(mat_imputed) <- dimnames(mat_norm)

cat("Missing after imputation:", sum(is.na(mat_imputed)), "\n")
# ------------------------------------------------------------------------------
# 6) QC plot: PCA after normalization + imputation (shown in R, not saved)
# ------------------------------------------------------------------------------
# PCA is performed on samples (columns). No scaling is applied because the data
# are already log-transformed and normalized.
pca_res <- prcomp(t(mat_imputed), center = TRUE, scale. = FALSE)

pca_df <- as_tibble(pca_res$x, rownames = "sample_id") %>%
  left_join(meta_lip, by = "sample_id")

var_expl <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)

p_pca <- ggplot(pca_df, aes(
  x = PC1, y = PC2, color = timepoint,
  shape = factor(any_brain_lesion)
)) +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_viridis_d(option = "D", end = 0.9) +
  theme_bw() +
  labs(
    title = "PCA – Lipidomics (post normalization + imputation)",
    x = paste0("PC1 (", var_expl[1], "%)"),
    y = paste0("PC2 (", var_expl[2], "%)")
  )

print(p_pca)

# ------------------------------------------------------------------------------
# 7) Save normalized + imputed matrix
# ------------------------------------------------------------------------------
expr_out <- bind_cols(
  expr_clean %>% select(all_of(info_cols)),
  as_tibble(mat_imputed) %>% setNames(bio_cols)
)

write_delim(expr_out, out_file, delim = "\t")

cat("\n------------------------------------------------\n")
cat("Saved normalized + imputed matrix to:", out_file, "\n")
cat("------------------------------------------------\n")
