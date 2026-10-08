# ==============================================================================
# Project: Preterm EV Multi-Omics – Proteomics
# Script:  05_prot_clustering_time_trajectories.R
# Author:  Nicole Grinovero
# Date:    Jan 2026
#
# Aim:
#   Build protein temporal trajectories (T0–T4), summarise each protein by the
#   median abundance within each timepoint, z-score each trajectory across time,
#   and perform k-means clustering.
#
# Optional QC (displayed in R, not saved):
#   - Elbow curve (total within-cluster SS) across a grid of k
#   - Average silhouette width across a grid of k
#
# Rationale for k:
#   - Elbow/silhouette diagnostics support a parsimonious macro-structure around
#     k = 2, corresponding to decreasing vs increasing temporal trajectories.
#   - For biological interpretability, we further subdivide these two macro
#     patterns into early vs late responders, yielding four subclusters (2×2).
#
# Inputs:
#   - data/proteomics/processed/proteins_norm_imputed.txt
#   - data/metadata/sample_metadata.txt
#   - outputs/tables/preterm/proteomics/proteins_LMM_time_results.txt
#
# Outputs:
#   - outputs/tables/preterm/proteomics/proteins_time_kmeans_clusters.txt
#   - outputs/tables/preterm/proteomics/prot_clustering_inputs.rds
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(tidyr)
})

# ------------------------------------------------------------------------------
# 0) Configuration
# ------------------------------------------------------------------------------
FDR_CUTOFF <- 0.05

SEED <- 12345
NSTART <- 50

# Final clustering used downstream
K_FINAL <- 4
cl_levels <- c("C1A", "C1B", "C2A", "C2B")

tp_levels <- c("T0", "T1", "T2", "T3", "T4")

# Optional QC: show elbow + silhouette in R (NOT saved)
DO_K_DIAGNOSTICS <- TRUE
K_GRID <- 2:10
NSTART_DIAG <- 25

expr_path <- here("data", "proteomics", "processed", "proteins_norm_imputed.txt")
meta_path <- here("data", "metadata", "sample_metadata.txt")
stats_path <- here("outputs", "tables", "preterm", "proteomics", "proteins_LMM_time_results.txt")

out_dir <- here("outputs", "tables", "preterm", "proteomics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_tbl <- file.path(out_dir, "proteins_time_kmeans_clusters.txt")
out_rds <- file.path(out_dir, "prot_clustering_inputs.rds")

# ------------------------------------------------------------------------------
# 1) Load
# ------------------------------------------------------------------------------
expr <- read_delim(expr_path, delim = "\t", show_col_types = FALSE)
meta <- read_delim(meta_path, delim = "\t", show_col_types = FALSE) %>%
  mutate(timepoint = factor(timepoint, levels = tp_levels))
stats <- read_delim(stats_path, delim = "\t", show_col_types = FALSE)

stopifnot(all(c("UniProtID", "GeneName") %in% names(expr)))
stopifnot(all(c("sample_id", "timepoint") %in% names(meta)))
stopifnot(all(c("UniProtID", "FDR_global") %in% names(stats)))

# Expression sample columns in metadata order 
bio_cols <- meta$sample_id[meta$sample_id %in% names(expr)]
if (length(bio_cols) == 0) stop("No sample columns matched between expr and metadata.")


# ------------------------------------------------------------------------------
# 2) Proteins with significant global time effect (LMM)
# ------------------------------------------------------------------------------
sig_ids <- stats %>%
  filter(FDR_global < FDR_CUTOFF) %>%
  distinct(UniProtID) %>%
  pull(UniProtID)

cat("Proteins significant for global time effect:", length(sig_ids), "\n")
if (length(sig_ids) < max(K_GRID)) {
  warning("Significant proteins (n=", length(sig_ids), ") may be small for larger k in K_GRID.")
}
if (length(sig_ids) < K_FINAL) {
  stop("Too few significant proteins for k-means clustering with K_FINAL=", K_FINAL)
}

# ------------------------------------------------------------------------------
# 3) Median trajectory matrix (T0-T4) and z-score across time
# ------------------------------------------------------------------------------
# Summarize repeated measures by median within each timepoint, per protein.
median_tp <- expr %>%
  select(UniProtID, GeneName, all_of(bio_cols)) %>%
  pivot_longer(all_of(bio_cols), names_to = "sample_id", values_to = "intensity") %>%
  left_join(meta %>% select(sample_id, timepoint), by = "sample_id") %>%
  filter(!is.na(timepoint), UniProtID %in% sig_ids) %>%
  group_by(UniProtID, GeneName, timepoint) %>%
  summarise(med = median(intensity, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = timepoint, values_from = med) %>%
  arrange(UniProtID)

mat_raw <- median_tp %>%
  select(all_of(tp_levels)) %>%
  as.matrix()
rownames(mat_raw) <- median_tp$UniProtID
colnames(mat_raw) <- tp_levels

# Z-score each protein across time (row-wise).
mat_z <- t(scale(t(mat_raw)))

# ------------------------------------------------------------------------------
# 4) OPTIONAL QC: elbow + silhouette (shown in R, not saved)
# ------------------------------------------------------------------------------
if (isTRUE(DO_K_DIAGNOSTICS)) {
  suppressPackageStartupMessages({
    library(cluster)  # silhouette
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
# 5) k-means clustering (K_FINAL) 
# ------------------------------------------------------------------------------
set.seed(SEED)
km <- kmeans(mat_z, centers = K_FINAL, nstart = NSTART, iter.max = 100)

# Centroid-derived features used to assign interpretable labels:
#   - family: decreasing vs increasing (based on T4 - T0 sign)
#   - late_score: magnitude of late change relative to early change
cent <- as.data.frame(km$centers) %>%
  tibble::rownames_to_column("cluster") %>%
  mutate(cluster = as.integer(cluster)) %>%
  transmute(
    cluster,
    delta_total = T4 - T0,
    delta_early = T1 - T0,
    delta_late  = T4 - T2,
    family      = if_else(delta_total < 0, "dec", "inc"),
    late_score  = abs(delta_late) - abs(delta_early) 
  )

n_dec <- sum(cent$family == "dec")
n_inc <- sum(cent$family == "inc")

# The labels below assume two decreasing and two increasing centroids (2 x 2)
if (n_dec != 2 || n_inc != 2) {
  stop("Expected 2 decreasing and 2 increasing clusters, got ", n_dec, " and ", n_inc, ".")
}


# Within each family, assign "late" first (B) then "early" (A)
dec_ids <- cent %>% filter(family == "dec") %>% arrange(desc(late_score)) %>% pull(cluster)
inc_ids <- cent %>% filter(family == "inc") %>% arrange(desc(late_score)) %>% pull(cluster)

map_tbl <- tibble(
  cluster = c(dec_ids, inc_ids),
  Cluster = c("C1B", "C1A", "C2B", "C2A")[seq_len(length(c(dec_ids, inc_ids)))]
)

cl <- map_tbl$Cluster[match(km$cluster, map_tbl$cluster)]
cl <- factor(cl, levels = cl_levels)

# ------------------------------------------------------------------------------
# 6) Row order for heatmap: Cluster then distance-to-centroid
# ------------------------------------------------------------------------------
centers_mat <- km$centers[km$cluster, , drop = FALSE]
dist_to_center <- rowSums((mat_z - centers_mat) ^ 2)

ord <- order(cl, dist_to_center)
mat_plot  <- mat_z[ord, , drop = FALSE]
row_split <- factor(as.character(cl[ord]), levels = cl_levels)

# ------------------------------------------------------------------------------
# 7) T4 vs T0 direction annotation
# ------------------------------------------------------------------------------
dir_tbl <- stats %>%
  distinct(UniProtID, .keep_all = TRUE) %>%
  transmute(
    UniProtID,
    T4vsT0 = case_when(
      !is.na(FDR_T4_T0) & FDR_T4_T0 < FDR_CUTOFF & beta_T4_vs_T0 > 0 ~ "Up",
      !is.na(FDR_T4_T0) & FDR_T4_T0 < FDR_CUTOFF & beta_T4_vs_T0 < 0 ~ "Down",
      TRUE ~ "NS"
    )
  )

ann_vec <- dir_tbl$T4vsT0[match(rownames(mat_plot), dir_tbl$UniProtID)]
ann_vec[is.na(ann_vec)] <- "NS"
ann_vec <- factor(ann_vec, levels = c("Down", "NS", "Up"))

# ------------------------------------------------------------------------------
# 8) Export cluster table
# ------------------------------------------------------------------------------
cluster_tbl <- tibble(
  UniProtID = rownames(mat_z),
  Cluster   = as.character(cl)
) %>%
  left_join(median_tp %>% select(UniProtID, GeneName), by = "UniProtID") %>%
  left_join(
    stats %>% select(UniProtID, PG.ProteinDescriptions, FDR_global, beta_T4_vs_T0, FDR_T4_T0),
    by = "UniProtID"
  ) %>%
  mutate(Cluster = factor(Cluster, levels = cl_levels)) %>%
  select(UniProtID, GeneName, PG.ProteinDescriptions, Cluster, FDR_global, beta_T4_vs_T0, FDR_T4_T0) %>%
  arrange(Cluster, UniProtID)

write_delim(cluster_tbl, out_tbl, delim = "\t")

# ------------------------------------------------------------------------------
# 9) Minimal bundle for 05_figures/fig4_proteome_remodeling.R
# ------------------------------------------------------------------------------
saveRDS(
  list(
    tp_levels = tp_levels,
    mat_plot  = mat_plot,
    row_split = row_split,
    ann_vec   = ann_vec
  ),
  out_rds
)

cat("\n[OK] Saved:\n - ", out_tbl, "\n - ", out_rds, "\n", sep = "")
