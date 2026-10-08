# ==============================================================================
# Project: Preterm EV Multi-Omics
# Figure:  Figure 4 – Proteomics developmental remodeling in plasma EVs
# Script:  fig4_proteome_remodeling.R
# Author:  Nicole Grinovero
# Date:    Jan 2026
#
# Panels:
#   A) Circular heatmap (z-scored trajectories T0–T4) + outer ring (T4 vs T0)
#   B) Cluster profiles plots (median ± IQR)
#   C) Volcano plot (T4 vs T0)
#   D) Representative proteins 
#   E) Enrichment analysis (GO:BP + Reactome; Up and Down)
#
# Inputs:
#   - outputs/tables/preterm/proteomics/prot_clustering_inputs.rds
#       list(tp_levels, mat_plot, row_split, ann_vec)
#   - outputs/tables/preterm/proteomics/proteins_LMM_time_results.txt
#   - data/proteomics/processed/proteins_norm_imputed.txt
#   - data/metadata/sample_metadata.txt
#   - outputs/tables/preterm/proteomics/enrich_GOBP_Reactome_T4vsT0.txt
#
# Outputs:
#   - outputs/figures/fig4_proteomics_remodeling/fig4a_circular_heatmap.svg
#   - outputs/figures/fig4_proteomics_remodeling/fig4b_cluster_profiles.svg
#   - outputs/figures/fig4_proteomics_remodeling/fig4c_volcano_T4_vs_T0.svg
#   - outputs/figures/fig4_proteomics_remodeling/fig4d_boxplots.svg
#   - outputs/figures/fig4_proteomics_remodeling/fig4e_enrichment_T4_vs_T0.svg
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(patchwork)
  library(svglite)
  
  # Panel A
  library(circlize)
  library(ComplexHeatmap)
  library(grid)
})

set.seed(1)

# ------------------------------------------------------------------------------
# 0) Global style
# ------------------------------------------------------------------------------
BASE_SIZE <- 13

theme_qomics <- function() {
  theme_bw(base_size = BASE_SIZE) +
    theme(
      text = element_text(color = "grey10"),
      axis.text = element_text(color = "grey10"),
      axis.title = element_text(color = "grey10"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 8, r = 6, b = 6, l = 6)
    )
}

save_svg <- function(plot, filename, width, height) {
  ggsave(
    filename = filename,
    plot     = plot,
    width    = width,
    height   = height,
    units    = "in",
    device   = svglite::svglite
  )
}

# ------------------------------------------------------------------------------
# 1) Paths
# ------------------------------------------------------------------------------
out_dir <- here("outputs", "figures", "fig4_proteomics_remodeling")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rds_cluster <- here("outputs","tables","preterm","proteomics","prot_clustering_inputs.rds")
stats_path  <- here("outputs","tables","preterm","proteomics","proteins_LMM_time_results.txt")
expr_path   <- here("data","proteomics","processed","proteins_norm_imputed.txt")
meta_path   <- here("data","metadata","sample_metadata.txt")
enrich_path <- here("outputs","tables","preterm","proteomics","enrich_GOBP_Reactome_T4vsT0.txt")

out_A <- file.path(out_dir, "fig4a_circular_heatmap.svg")
out_B <- file.path(out_dir, "fig4b_cluster_profiles.svg")
out_C <- file.path(out_dir, "fig4c_volcano_T4_vs_T0.svg")
out_D <- file.path(out_dir, "fig4d_boxplots.svg")
out_E <- file.path(out_dir, "fig4e_enrichment_T4_vs_T0.svg")

stopifnot(
  file.exists(rds_cluster),
  file.exists(stats_path),
  file.exists(expr_path),
  file.exists(meta_path),
  file.exists(enrich_path)
)

# ------------------------------------------------------------------------------
# 2) Global thresholds + palette (used consistently across panels)
# ------------------------------------------------------------------------------
FDR_CUTOFF <- 0.05
FC_CUTOFF  <- log2(1.5)

COL_UP   <- "#7F4E57"    # burgundy
COL_DOWN <- "#00344AFF"  # teal/blue
COL_NS   <- "grey75"
COL_NS2  <- "#D0D0D0"

# Z-score color scale (Panels A/B)
CLIP_Z <- 2
pal_z <- circlize::colorRamp2(
  c(-CLIP_Z, 0, CLIP_Z),
  c("#3B6FB6", "#F7F7F7", "#B33A3A")
)

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------
blend_white <- function(col, p = 0.5) {
  grDevices::colorRampPalette(c(col, "#FFFFFF"))(101)[round(p * 100) + 1]
}

plot_lollipop <- function(df, fill_col) {
  df <- df %>%
    arrange(neglog10) %>%
    mutate(
      Term = str_wrap(Term, width = 32),
      Term = factor(Term, levels = Term)
    )
  
  ggplot(df, aes(x = neglog10, y = Term)) +
    geom_segment(aes(x = 0, xend = neglog10, yend = Term),
                 linewidth = 1.1, color = "black") +
    geom_point(aes(size = Count),
               shape = 21, color = "grey70", fill = fill_col,
               alpha = 0.9, stroke = 1.1) +
    scale_size_continuous(range = c(2, 5), name = "Protein count") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.10))) +
    labs(x = expression(paste("-log"[10],"(FDR)")), y = NULL) +
    theme_qomics() +
    theme(
      panel.grid.major.x = element_line(linetype = "dashed", linewidth = 0.3, color = "grey80"),
      panel.grid.major.y = element_line(linewidth = 0.25, color = "grey92"),
      panel.grid.minor   = element_blank(),
      strip.background   = element_rect(fill = "grey95", color = "black", linewidth = 0.25),
      strip.text         = element_text(face = "bold", size = BASE_SIZE - 2),
      legend.position    = "right"
    )
}

# ------------------------------------------------------------------------------
# 3) Load clustering bundle (produced upstream)
# ------------------------------------------------------------------------------
clu <- readRDS(rds_cluster)
stopifnot(all(c("tp_levels","mat_plot","row_split","ann_vec") %in% names(clu)))

tp_levels <- clu$tp_levels
mat_z_ord <- clu$mat_plot
row_split <- clu$row_split
ann_vec   <- clu$ann_vec

stopifnot(
  all(tp_levels %in% colnames(mat_z_ord)),
  nrow(mat_z_ord) == length(row_split),
  nrow(mat_z_ord) == length(ann_vec)
)

mat_z_clip <- pmax(pmin(mat_z_ord, CLIP_Z), -CLIP_Z)

# ==============================================================================
# PANEL A – Circular heatmap + outer ring annotation (T4 vs T0)
# ==============================================================================
mat_circ  <- mat_z_clip[, rev(tp_levels), drop = FALSE]
split_fac <- factor(as.character(row_split), levels = levels(row_split))

lev0 <- levels(split_fac)
if (!any(grepl("^Cluster", lev0))) {
  levels(split_fac) <- gsub("^C", "Cluster ", lev0)
}

t4_cols <- c("Down" = COL_DOWN, "NS" = COL_NS2, "Up" = COL_UP)
ann_fix <- as.character(ann_vec)
ann_fix[is.na(ann_fix)] <- "NS"
ann_fix <- factor(ann_fix, levels = c("Down","NS","Up"))

svglite::svglite(out_A, width = 6.2, height = 6.2)

circlize::circos.clear()
circlize::circos.par(
  start.degree = 90,
  gap.after = c(rep(2, length(levels(split_fac)) - 1), 12),
  cell.padding = c(0, 0, 0, 0),
  track.margin = c(0.003, 0.010),
  points.overflow.warning = FALSE
)

circlize::circos.heatmap(
  mat_circ,
  split = split_fac,
  col = pal_z,
  cluster = TRUE,
  dend.side = "inside",
  dend.track.height = 0.15,
  show.sector.labels = TRUE,
  track.height = 0.34
)

circlize::circos.trackPlotRegion(
  track.index = 1,
  bg.border = NA,
  panel.fun = function(x, y) {
    circlize::circos.rect(
      xleft   = CELL_META$cell.xlim[1],
      xright  = CELL_META$cell.xlim[2],
      ybottom = CELL_META$cell.ylim[1],
      ytop    = CELL_META$cell.ylim[2],
      col     = NA,
      border  = "black",
      lwd     = 0.8
    )
  }
)

first_sector <- levels(split_fac)[1]
circlize::circos.trackPlotRegion(
  track.index = get.current.track.index(),
  bg.border = NA,
  panel.fun = function(x, y) {
    if (CELL_META$sector.index == first_sector) {
      n_tp  <- ncol(mat_circ)
      atpos <- seq_len(n_tp) - 0.5
      circos.yaxis(
        side = "left",
        at = atpos,
        labels = rev(colnames(mat_circ)),
        labels.cex = 0.9,
        tick.length = 0.02,
        labels.niceFacing = TRUE
      )
    }
  }
)

circlize::circos.trackPlotRegion(
  ylim = c(0, 1),
  bg.border = NA,
  track.height = 0.055,
  panel.fun = function(x, y) {
    idx_global    <- CELL_META$subset
    ord_in_sector <- CELL_META$row_order
    idx_plot      <- idx_global[ord_in_sector]
    
    cats <- ann_fix[idx_plot]
    n <- length(idx_plot)
    
    circlize::circos.rect(
      xleft   = seq_len(n) - 1,
      xright  = seq_len(n),
      ybottom = 0,
      ytop    = 1,
      col     = unname(t4_cols[as.character(cats)]),
      border  = NA
    )
  }
)

lgd_z <- ComplexHeatmap::Legend(
  title = "Z-score",
  col_fun = pal_z,
  at = c(-CLIP_Z, 0, CLIP_Z),
  direction = "horizontal"
)
lgd_t4 <- ComplexHeatmap::Legend(
  title = "T4 vs T0",
  labels = c("Down at TEA", "Not significant", "Up at TEA"),
  legend_gp = grid::gpar(fill = unname(t4_cols[c("Down","NS","Up")])),
  direction = "horizontal"
)

ComplexHeatmap::draw(
  ComplexHeatmap::packLegend(
    lgd_z, lgd_t4,
    direction = "horizontal",
    gap = grid::unit(6, "mm")
  ),
  x = grid::unit(0.5, "npc"),
  y = grid::unit(0.5, "npc"),
  just = c("center", "center")
)

dev.off()
circlize::circos.clear()

# ==============================================================================
# PANEL B – Cluster profile plots 
# ==============================================================================
COL_1B <- COL_DOWN
COL_2B <- COL_UP
COL_1A <- blend_white(COL_1B, p = 0.5)
COL_2A <- blend_white(COL_2B, p = 0.5)

cluster_cols <- c(
  "Cluster 1A" = COL_1A, "Cluster 1B" = COL_1B,
  "Cluster 2A" = COL_2A, "Cluster 2B" = COL_2B
)

trend_df <- mat_z_ord %>%
  as.data.frame() %>%
  tibble::rownames_to_column("UniProtID") %>%
  mutate(ClusterRaw = as.character(row_split)) %>%
  pivot_longer(cols = all_of(tp_levels), names_to = "Timepoint", values_to = "z") %>%
  mutate(Timepoint = factor(Timepoint, levels = tp_levels))

trend_sum <- trend_df %>%
  group_by(ClusterRaw, Timepoint) %>%
  summarise(
    med = median(z, na.rm = TRUE),
    q1  = quantile(z, 0.25, na.rm = TRUE),
    q3  = quantile(z, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

cl_delta <- trend_sum %>%
  filter(Timepoint %in% c("T0","T4")) %>%
  select(ClusterRaw, Timepoint, med) %>%
  pivot_wider(names_from = Timepoint, values_from = med) %>%
  mutate(delta = T4 - T0)

dec_raw <- cl_delta %>% arrange(delta) %>% slice(1:2) %>% pull(ClusterRaw)
inc_raw <- cl_delta %>% arrange(desc(delta)) %>% slice(1:2) %>% pull(ClusterRaw)

dec_tbl <- cl_delta %>% filter(ClusterRaw %in% dec_raw) %>% mutate(abs_delta = abs(delta)) %>% arrange(abs_delta)
inc_tbl <- cl_delta %>% filter(ClusterRaw %in% inc_raw) %>% mutate(abs_delta = abs(delta)) %>% arrange(abs_delta)

map_labels <- c("Cluster 1A","Cluster 1B","Cluster 2A","Cluster 2B")
names(map_labels) <- c(
  as.character(dec_tbl$ClusterRaw[1]),
  as.character(dec_tbl$ClusterRaw[2]),
  as.character(inc_tbl$ClusterRaw[1]),
  as.character(inc_tbl$ClusterRaw[2])
)

trend_df <- trend_df %>%
  mutate(
    Cluster = unname(map_labels[as.character(ClusterRaw)]),
    Cluster = factor(Cluster, levels = c("Cluster 1A","Cluster 1B","Cluster 2A","Cluster 2B"))
  ) %>% filter(!is.na(Cluster))

trend_sum <- trend_sum %>%
  mutate(
    Cluster = unname(map_labels[as.character(ClusterRaw)]),
    Cluster = factor(Cluster, levels = c("Cluster 1A","Cluster 1B","Cluster 2A","Cluster 2B"))
  ) %>% filter(!is.na(Cluster))

p_B <- ggplot() +
  geom_line(
    data = trend_df,
    aes(x = Timepoint, y = z, group = UniProtID, color = Cluster),
    alpha = 0.09, linewidth = 0.32
  ) +
  geom_ribbon(
    data = trend_sum,
    aes(x = Timepoint, ymin = q1, ymax = q3, fill = Cluster),
    alpha = 0.18
  ) +
  geom_line(
    data = trend_sum,
    aes(x = Timepoint, y = med, group = Cluster, color = Cluster),
    linewidth = 2.2
  ) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.5, alpha = 0.5) +
  facet_wrap(~ Cluster, ncol = 2, scales = "fixed") +
  coord_cartesian(ylim = c(-2, 2)) +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) +
  scale_fill_manual(values = cluster_cols) +
  scale_color_manual(values = cluster_cols) +
  labs(x = NULL, y = "Protein abundance (z-score)") +
  theme_qomics() +
  theme(
    panel.grid = element_blank(),
    panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.background = element_blank(),
    strip.text       = element_text(size = BASE_SIZE),
    legend.position  = "none"
  )

print(p_B)

save_svg(p_B, out_B, width = 4.8, height = 4.0)

# ==============================================================================
# PANEL C – Volcano plot (T4 vs T0)
# ==============================================================================
stats <- read_delim(stats_path, delim = "\t", show_col_types = FALSE)
stopifnot(all(c("GeneName","beta_T4_vs_T0","FDR_T4_T0") %in% names(stats)))

label_genes <- c("AFP","LDLR","CPVL","LMAN2","PGRMC1","MTTP","HBD","TOR1A","SHBG","APOA4","C9")

volcano_df <- stats %>%
  mutate(
    GeneName    = str_trim(GeneName),
    FDR_T4_T0   = pmax(FDR_T4_T0, 1e-300),
    neglog10_fdr = -log10(FDR_T4_T0),
    endpoint_class = case_when(
      FDR_T4_T0 < FDR_CUTOFF & beta_T4_vs_T0 >  FC_CUTOFF ~ "Up",
      FDR_T4_T0 < FDR_CUTOFF & beta_T4_vs_T0 < -FC_CUTOFF ~ "Down",
      TRUE ~ "NS"
    ),
    endpoint_class = factor(endpoint_class, levels = c("Down","NS","Up")),
    label = if_else(GeneName %in% label_genes, GeneName, NA_character_)
  )

cap_val <- as.numeric(quantile(volcano_df$neglog10_fdr, 0.995, na.rm = TRUE))
volcano_df <- volcano_df %>%
  mutate(
    sig_cap  = pmin(neglog10_fdr, cap_val),
    pt_size  = rescale(sig_cap, to = c(1.5, 3.5), from = range(sig_cap, na.rm = TRUE)),
    pt_alpha = rescale(sig_cap, to = c(0.2, 0.9), from = range(sig_cap, na.rm = TRUE))
  )

up_count   <- sum(volcano_df$endpoint_class == "Up",   na.rm = TRUE)
down_count <- sum(volcano_df$endpoint_class == "Down", na.rm = TRUE)

y_thr <- -log10(FDR_CUTOFF)
max_y <- max(volcano_df$neglog10_fdr, na.rm = TRUE)
min_x <- min(volcano_df$beta_T4_vs_T0, na.rm = TRUE)
max_x <- max(volcano_df$beta_T4_vs_T0, na.rm = TRUE)

p_C <- ggplot(volcano_df, aes(x = beta_T4_vs_T0, y = neglog10_fdr)) +
  geom_hline(yintercept = y_thr, linetype = "dashed",
             color = alpha("grey55", 0.65), linewidth = 0.4) +
  geom_vline(xintercept = c(-FC_CUTOFF, FC_CUTOFF), linetype = "dashed",
             color = alpha("grey55", 0.65), linewidth = 0.4) +
  geom_point(
    data = volcano_df %>% filter(endpoint_class == "NS"),
    aes(size = pt_size, alpha = pt_alpha),
    color = COL_NS
  ) +
  geom_point(
    data = volcano_df %>% filter(endpoint_class != "NS"),
    aes(color = endpoint_class, size = pt_size, alpha = pt_alpha)
  ) +
  scale_color_manual(values = c("Up" = COL_UP, "Down" = COL_DOWN)) +
  scale_size_identity() +
  scale_alpha_identity() +
  geom_text_repel(
    data = volcano_df %>% filter(!is.na(label)),
    aes(label = label),
    box.padding = 0.55,
    point.padding = 0.25,
    size = 3.3,
    max.overlaps = Inf,
    segment.color = "grey70",
    min.segment.length = 0
  ) +
  annotate("text", x = min_x, y = -log10(FDR_CUTOFF) + 0.9, label = paste0("▼ ", down_count),
           color = COL_DOWN, fontface = "bold", size = 5.0, hjust = 0) +
  annotate("text", x = max_x, y = -log10(FDR_CUTOFF) + 0.9, label = paste0("▲ ", up_count),
           color = COL_UP, fontface = "bold", size = 5.0, hjust = 1) +
  labs(
    x = "Effect size (T4 vs T0)",
    y = expression("-log"[10]*"(FDR)")
  ) +
  coord_cartesian(ylim = c(0, max_y + 0.9)) +
  theme_qomics() +
  theme(panel.grid = element_blank(), legend.position = "none")

print(p_C)

save_svg(p_C, out_C, width = 4.0, height = 4.0)

# ==============================================================================
# PANEL D – Representative proteins (boxplots)
# ==============================================================================
genes_keep <- c("AFP","RPS6","EIF4A2","MTTP","C9","PROS1")
tp5 <- c("T0","T1","T2","T3","T4")

expr <- read_delim(expr_path, delim = "\t", show_col_types = FALSE) %>%
  rename(Gene = GeneName) %>%
  mutate(Gene = str_trim(Gene))

meta <- read_delim(meta_path, delim = "\t", show_col_types = FALSE) %>%
  mutate(
    sample_id  = as.character(sample_id),
    patient_id = as.character(patient_id),
    timepoint  = factor(as.character(timepoint), levels = tp5)
  )

sample_cols <- intersect(colnames(expr), unique(meta$sample_id))
stopifnot(length(sample_cols) > 0)

df_D <- expr %>%
  select(Gene, all_of(sample_cols)) %>%
  filter(Gene %in% genes_keep) %>%
  pivot_longer(cols = -Gene, names_to = "sample_id", values_to = "intensity") %>%
  left_join(meta, by = "sample_id") %>%
  filter(!is.na(patient_id), !is.na(timepoint), !is.na(intensity)) %>%
  mutate(
    Direction = case_when(
      Gene %in% c("MTTP","C9","PROS1") ~ "Up",
      Gene %in% c("AFP","RPS6","EIF4A2") ~ "Down",
      TRUE ~ NA_character_
    ),
    Direction = factor(Direction, levels = c("Down","Up")),
    Gene = factor(Gene, levels = genes_keep)
  )

p_D <- ggplot(df_D, aes(x = timepoint, y = intensity)) +
  geom_jitter(width = 0.08, size = 1.4, alpha = 0.6, color = "grey70") +
  geom_boxplot(aes(fill = Direction), width = 0.55, outlier.shape = NA,
               alpha = 0.6, linewidth = 0.6, color = "grey35") +
  stat_summary(fun = median, geom = "point", shape = 21, size = 1.8,
               stroke = 0.3, fill = "white", color = "grey35") +
  scale_fill_manual(values = c("Down" = COL_DOWN, "Up" = COL_UP)) +
  facet_wrap(~ Gene, ncol = 3, scales = "free_y") +
  labs(x = NULL, y = expression(paste("log"[2],"(abundance)")), fill = NULL) +
  theme_qomics() +
  theme(
    strip.text = element_text(size = 12),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.spacing = unit(0.9, "lines")
  )

print(p_D)

save_svg(p_D, out_D, width = 7.0, height = 4)

# ==============================================================================
# PANEL E – GO:BP + Reactome enrichment analysis
# ==============================================================================
enr <- read_delim(enrich_path, delim = "\t", show_col_types = FALSE)
stopifnot(all(c("Direction","Source","Term","FDR","Count","neglog10") %in% names(enr)))

df_plot <- enr %>%
  mutate(
    Direction = factor(Direction, levels = c("Down","Up")),
    Source    = factor(Source, levels = c("GO:BP","Reactome")),
    Term      = str_replace_all(Term, "_", " ")
  )

# -----------------------------
# Downregulated proteins
# -----------------------------
df_down <- df_plot %>%
  filter(Direction == "Down") %>%
  arrange(neglog10) %>%
  mutate(
    Term = str_wrap(Term, width = 36),
    Term = factor(Term, levels = Term)
  )

p_E_down <- ggplot(df_down, aes(x = neglog10, y = Term)) +
  geom_point(aes(size = Count), color = COL_DOWN, alpha = 0.7) +
  scale_size(range = c(2, 6.5), name = "Protein count") +
  scale_x_continuous(
    expand = expansion(mult = c(0.2, 0.2)),
    name = expression(paste("-log"[10], "(FDR)"))
  ) +
  labs(y = NULL) +
  theme_qomics() +
  theme(
    axis.text.y = element_text(size = BASE_SIZE - 1),
    panel.grid.major.x = element_line(linetype = "dashed", linewidth = 0.3, color = "grey80"),
    panel.grid.major.y = element_line(linewidth = 0.25, color = "grey92"),
    panel.grid.minor   = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = BASE_SIZE)
  )

# -----------------------------
# Upregulated proteins
# -----------------------------
df_up <- df_plot %>%
  filter(Direction == "Up") %>%
  arrange(neglog10) %>%
  mutate(
    Term = str_wrap(Term, width = 36),
    Term = factor(Term, levels = Term)
  )

p_E_up <- ggplot(df_up, aes(x = neglog10, y = Term)) +
  geom_point(aes(size = Count), color = COL_UP, alpha = 0.7) +
  scale_size(range = c(2, 6.5), name = "Protein count") +
  scale_x_continuous(
    expand = expansion(mult = c(0.2, 0.2)),
    name = expression(paste("-log"[10], "(FDR)"))
  ) +
  labs(y = NULL) +
  theme_qomics() +
  theme(
    axis.text.y = element_text(size = BASE_SIZE - 1),
    panel.grid.major.x = element_line(linetype = "dashed", linewidth = 0.3, color = "grey80"),
    panel.grid.major.y = element_line(linewidth = 0.25, color = "grey92"),
    panel.grid.minor   = element_blank(),
    legend.position = "right",   
    plot.title = element_text(face = "bold", size = BASE_SIZE)
  )

# Combine
p_E <- p_E_down + p_E_up + plot_layout(ncol = 2)

print(p_E)

save_svg(p_E, out_E, width = 11.5, height = 2.6)

# ------------------------------------------------------------------------------
# Done
# ------------------------------------------------------------------------------
cat(
  "\n[OK] Figure 4 panels saved in:\n",
  " - ", out_A, "\n",
  " - ", out_B, "\n",
  " - ", out_C, "\n",
  " - ", out_D, "\n",
  " - ", out_E, "\n",
  sep = ""
)
