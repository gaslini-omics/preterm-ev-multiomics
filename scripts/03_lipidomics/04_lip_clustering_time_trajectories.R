# ==============================================================================
# Project: Preterm EV Multi-Omics – Lipidomics
# Script:  04_lip_clustering_time_trajectories.R
# Author:  Nicole Grinovero
# Date:    Jan 2026
#
# Aim:
#   Build lipid temporal trajectories (T0–T4), summarise each lipid by the
#   median abundance within each timepoint, z-score each trajectory across time,
#   and perform k-means clustering (K_FINAL = 4).
#
# Optional QC (displayed in R, not saved):
#   - Elbow curve (total within-cluster SS) across a grid of k
#   - Average silhouette width across a grid of k
#
# Inputs:
#   - data/lipidomics/processed/lipids_norm_imputed.txt
#   - data/metadata/sample_metadata.txt
#   - outputs/tables/preterm/lipidomics/lipids_LMM_time_results.txt
#
# Outputs:
#   - outputs/tables/preterm/lipidomics/lipids_time_kmeans_clusters.txt
#   - outputs/tables/preterm/lipidomics/lip_clustering_inputs.rds
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# ------------------------------------------------------------------------------
# 0) Configuration
# ------------------------------------------------------------------------------
FDR_CUTOFF <- 0.05

SEED    <- 12345
NSTART  <- 50
K_FINAL <- 4

tp_levels <- c("T0", "T1", "T2", "T3", "T4")

# Row-wise z-score clipping
CLIP_Z <- 2

# Optional QC: show elbow + silhouette in R (NOT saved)
DO_K_DIAGNOSTICS <- TRUE
K_GRID <- 2:10
NSTART_DIAG <- 25

# Drop deuterated internal standards ("(dX)" in the name)
DROP_ISOTOPES <- TRUE

# Paths
expr_path  <- here("data", "lipidomics", "processed", "lipids_norm_imputed.txt")
meta_path  <- here("data", "metadata", "sample_metadata.txt")
stats_path <- here("outputs", "tables", "preterm", "lipidomics", "lipids_LMM_time_results.txt")

out_dir <- here("outputs", "tables", "preterm", "lipidomics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_tbl <- file.path(out_dir, "lipids_time_kmeans_clusters.txt")
out_rds <- file.path(out_dir, "lip_clustering_inputs.rds")

# ------------------------------------------------------------------------------
# 1) Load + basic checks
# ------------------------------------------------------------------------------
stopifnot(file.exists(expr_path), file.exists(meta_path), file.exists(stats_path))

expr <- read_delim(expr_path, delim = "\t", show_col_types = FALSE)

# Required columns in expression table
stopifnot("LipidName" %in% names(expr))
stopifnot("Class" %in% names(expr))

# Deuterated internal standards
if (isTRUE(DROP_ISOTOPES)) {
  expr <- expr %>% filter(!str_detect(LipidName, "\\(d[0-9]+\\)"))
}

meta <- read_delim(meta_path, delim = "\t", show_col_types = FALSE)
stopifnot(all(c("sample_id", "timepoint") %in% names(meta)))

meta <- meta %>%
  mutate(
    sample_id = as.character(sample_id),
    timepoint = factor(as.character(timepoint), levels = tp_levels)
  )

stats <- read_delim(stats_path, delim = "\t", show_col_types = FALSE)
stopifnot(all(c("LipidName", "FDR_global", "beta_T4_vs_T0", "FDR_T4_T0") %in% names(stats)))

# Expression sample columns in metadata order
bio_cols <- meta$sample_id[meta$sample_id %in% names(expr)]
if (length(bio_cols) == 0) stop("No sample columns matched between expr and metadata.")

# ------------------------------------------------------------------------------
# 2) Lipids with significant global time effect (LMM)
# ------------------------------------------------------------------------------
sig_ids <- stats %>%
  filter(FDR_global < FDR_CUTOFF) %>%
  distinct(LipidName) %>%
  pull(LipidName)

cat("Lipids significant for global time effect:", length(sig_ids), "\n")
if (length(sig_ids) == 0) stop("No lipids pass FDR_global cutoff. Nothing to cluster.")

if (isTRUE(DO_K_DIAGNOSTICS) && length(sig_ids) < max(K_GRID)) {
  warning("Significant lipids (n=", length(sig_ids), ") may be small for larger k in K_GRID.")
}

# ------------------------------------------------------------------------------
# 3) Median trajectory matrix (T0–T4) + row-wise z-score (and clipping)
# ------------------------------------------------------------------------------
# Keep only what we need from expr to avoid carrying extra columns
expr_long <- expr %>%
  select(LipidName, Class, all_of(bio_cols)) %>%
  pivot_longer(all_of(bio_cols), names_to = "sample_id", values_to = "intensity") %>%
  left_join(meta %>% select(sample_id, timepoint), by = "sample_id") %>%
  filter(!is.na(timepoint), LipidName %in% sig_ids)

# Median per lipid per timepoint
median_tp <- expr_long %>%
  group_by(LipidName, timepoint) %>%
  summarise(med = median(intensity, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = timepoint, values_from = med) %>%
  arrange(LipidName)

# Ensure all timepoints exist as columns (even if some are all NA)
missing_tp_cols <- setdiff(tp_levels, names(median_tp))
if (length(missing_tp_cols) > 0) {
  for (tt in missing_tp_cols) median_tp[[tt]] <- NA_real_
  median_tp <- median_tp %>% select(LipidName, all_of(tp_levels))
}

mat_raw <- median_tp %>% select(all_of(tp_levels)) %>% as.matrix()
rownames(mat_raw) <- median_tp$LipidName
colnames(mat_raw) <- tp_levels

# Row-wise z-score + clip
mat_z <- t(scale(t(mat_raw)))
mat_z <- pmax(pmin(mat_z, CLIP_Z), -CLIP_Z)

# ------------------------------------------------------------------------------
# 4) OPTIONAL QC: elbow + silhouette (shown in R, not saved)
# ------------------------------------------------------------------------------
if (isTRUE(DO_K_DIAGNOSTICS)) {
  suppressPackageStartupMessages({
    library(cluster)
    library(ggplot2)
  })
  
  set.seed(SEED)
  dmat <- dist(mat_z)
  
  diag_tbl <- lapply(K_GRID, function(k) {
    kmk <- kmeans(mat_z, centers = k, nstart = NSTART_DIAG, iter.max = 100)
    sil <- cluster::silhouette(kmk$cluster, dmat)
    tibble(
      k = k,
      tot_withinss = kmk$tot.withinss,
      avg_sil = mean(sil[, "sil_width"])
    )
  }) %>% bind_rows()
  
  cat("\n[QC] k-means diagnostics (not saved):\n")
  print(diag_tbl)
  
  p_elbow <- ggplot(diag_tbl, aes(k, tot_withinss)) +
    geom_line(linewidth = 0.7, color = "grey35") +
    geom_point(size = 2.0, color = "grey20") +
    scale_x_continuous(breaks = K_GRID) +
    labs(x = "k", y = "Total within-cluster SS", title = "Elbow method") +
    theme_bw(base_size = 11)
  
  p_sil <- ggplot(diag_tbl, aes(k, avg_sil)) +
    geom_line(linewidth = 0.7, color = "grey35") +
    geom_point(size = 2.0, color = "grey20") +
    scale_x_continuous(breaks = K_GRID) +
    labs(x = "k", y = "Average silhouette width", title = "Silhouette") +
    theme_bw(base_size = 11)
  
  print(p_elbow)
  print(p_sil)
}

# ------------------------------------------------------------------------------
# 5) k-means clustering (K_FINAL = 4)
# ------------------------------------------------------------------------------
set.seed(SEED)
km <- kmeans(mat_z, centers = K_FINAL, nstart = NSTART, iter.max = 100)

cl <- factor(as.integer(km$cluster), levels = seq_len(K_FINAL))

# ------------------------------------------------------------------------------
# 6) Row order for plotting: Cluster then distance-to-centroid
# ------------------------------------------------------------------------------
centers_mat <- km$centers[km$cluster, , drop = FALSE]
dist_to_center <- rowSums((mat_z - centers_mat) ^ 2)

ord <- order(cl, dist_to_center)

mat_plot  <- mat_z[ord, , drop = FALSE]
row_split <- factor(as.character(cl[ord]), levels = as.character(seq_len(K_FINAL)))

# ------------------------------------------------------------------------------
# 7) T4 vs T0 direction annotation (Up/Down/NS) aligned to mat_plot
# ------------------------------------------------------------------------------
dir_tbl <- stats %>%
  distinct(LipidName, .keep_all = TRUE) %>%
  transmute(
    LipidName,
    T4vsT0 = case_when(
      !is.na(FDR_T4_T0) & FDR_T4_T0 < FDR_CUTOFF & beta_T4_vs_T0 > 0 ~ "Up",
      !is.na(FDR_T4_T0) & FDR_T4_T0 < FDR_CUTOFF & beta_T4_vs_T0 < 0 ~ "Down",
      TRUE ~ "NS"
    )
  )

ann_vec <- dir_tbl$T4vsT0[match(rownames(mat_plot), dir_tbl$LipidName)]
ann_vec[is.na(ann_vec)] <- "NS"
ann_vec <- factor(ann_vec, levels = c("Down", "NS", "Up"))

# ------------------------------------------------------------------------------
# 8) Class vector aligned to mat_plot
# ------------------------------------------------------------------------------
class_map <- expr %>% distinct(LipidName, Class)

class_vec <- class_map$Class[match(rownames(mat_plot), class_map$LipidName)]
class_vec[is.na(class_vec) | class_vec == ""] <- "Other"
class_vec <- as.character(class_vec)

# ------------------------------------------------------------------------------
# 9) Export cluster table 
# ------------------------------------------------------------------------------
cluster_tbl <- tibble(
  LipidName = rownames(mat_z),
  Cluster   = factor(as.integer(cl), levels = seq_len(K_FINAL))
) %>%
  left_join(dir_tbl, by = "LipidName") %>%
  left_join(
    stats %>% select(LipidName, FDR_global, beta_T4_vs_T0, FDR_T4_T0),
    by = "LipidName"
  ) %>%
  left_join(class_map, by = "LipidName") %>%
  mutate(
    Class = if_else(is.na(Class) | Class == "", "Other", as.character(Class))
  ) %>%
  select(LipidName, Class, Cluster, T4vsT0, FDR_global, beta_T4_vs_T0, FDR_T4_T0) %>%
  arrange(Cluster, LipidName)

write_delim(cluster_tbl, out_tbl, delim = "\t")

# ------------------------------------------------------------------------------
# 10) Minimal bundle for 05_figures/fig5_lipidome_remodeling.R
#     (everything already aligned to mat_plot order)
# ------------------------------------------------------------------------------
saveRDS(
  list(
    tp_levels = tp_levels,
    CLIP_Z    = CLIP_Z,
    K_FINAL   = K_FINAL,
    mat_plot  = mat_plot,
    row_split = row_split,
    ann_vec   = ann_vec,
    class_vec = class_vec,
    centers   = km$centers
  ),
  out_rds
)

cat(
  "\n[OK] Saved:\n",
  " - ", out_tbl, "\n",
  " - ", out_rds, "\n",
  sep = ""
)

