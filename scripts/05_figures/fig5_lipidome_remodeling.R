# ==============================================================================
# Project: Preterm EV Multi-Omics – Lipidomics
# Figure:  Figure 5 – Lipidomics developmental remodeling in plasma EVs
# Script:  fig5_lipidome_remodeling.R
# Author:  Nicole Grinovero
# Date:    Jan 2026
#
# Panels:
#   A) Temporal heatmap: rows = timepoints (T0–T4), cols = lipids
#   B) Cluster profiles (median ± IQR) 
#   C) Volcano plot (T4 vs T0)
#   D) LION enrichment barplot 
#   E) Acyl-chain structural space for selected membrane lipid classes
#
# Inputs:
#   - outputs/tables/preterm/lipidomics/lip_clustering_inputs.rds
#   - outputs/tables/preterm/lipidomics/lipids_LMM_time_results.txt
#   - outputs/tables/preterm/lipidomics/LION-enrichment-report-T4-vs-T0/...
#
# Outputs:
#   - outputs/figures/fig5_lipidomics_remodeling/fig5a_lip_heatmap_k4.svg
#   - outputs/figures/fig5_lipidomics_remodeling/fig5b_lip_cluster_profiles.svg
#   - outputs/figures/fig5_lipidomics_remodeling/fig5c_lip_volcano_T4_vs_T0.svg
#   - outputs/figures/fig5_lipidomics_remodeling/fig5d_lip_LION_enrich_plot.svg
#   - outputs/figures/fig5_lipidomics_remodeling/fig5e_lip_acylchain_T4vsT0.svg
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(svglite)
  library(patchwork)
  
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  
  library(ggnewscale)
})

# ------------------------------------------------------------------------------
# 0) GLOBAL CONFIG
# ------------------------------------------------------------------------------

set.seed(1)

BASE_SIZE  <- 13
FDR_CUTOFF <- 0.05
FC_CUTOFF  <- log2(1.5)

# Palette / colors 
COL_UP   <- "#A15D5D"
COL_DOWN <- "#2F7F7B"
COL_NS   <- "grey75"
COL_NS2  <- "#D9D9D9"

# Heatmap z-score clipping + palette
CLIP_Z <- 2
pal_z <- circlize::colorRamp2(
  c(-CLIP_Z, 0, CLIP_Z),
  c("#3B6FB6", "#F7F7F7", "#B33A3A")
)

# Direction annotation colors (Panel A)
t4_cols <- c(
  `Down at TEA`     = COL_DOWN,
  `Not significant` = COL_NS2,
  `Up at TEA`       = COL_UP
)

# Heatmap fine-tuning
ANNO_THICK_MM   <- 2.2   # height of the top annotation
COLUMN_GAP_MM   <- 1.4
HEATMAP_BORDER  <- TRUE

# ------------------------------------------------------------------------------
# 1) PATHS
# ------------------------------------------------------------------------------

clust_path <- here("outputs", "tables", "preterm", "lipidomics", "lip_clustering_inputs.rds")
stats_path <- here("outputs", "tables", "preterm", "lipidomics", "lipids_LMM_time_results.txt")

lion_path  <- here(
  "outputs", "tables", "preterm", "lipidomics",
  "LION-enrichment-report-T4-vs-T0",
  "LION_enrichment_ranking_beta_T4_vs_T0.csv"
)

fig_dir <- here("outputs", "figures", "fig5_lipidomics_remodeling")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

out_A <- file.path(fig_dir, "fig5a_lip_heatmap_k4.svg")
out_B <- file.path(fig_dir, "fig5b_lip_cluster_profiles.svg")
out_C <- file.path(fig_dir, "fig5c_lip_volcano_T4_vs_T0.svg")
out_D <- file.path(fig_dir, "fig5d_lip_LION_enrich_plot.svg")
out_E <- file.path(fig_dir, "fig5e_lip_acylchain_T4vsT0.svg")

stopifnot(file.exists(clust_path), file.exists(stats_path))

# ------------------------------------------------------------------------------
# 2) HELPERS (STYLE + EXPORT)
# ------------------------------------------------------------------------------

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

save_ht_svg <- function(ht, filename, width, height, legends = NULL,
                        heatmap_legend_side = "right",
                        annotation_legend_side = "right") {
  svglite::svglite(filename, width = width, height = height)
  ComplexHeatmap::draw(
    ht,
    annotation_legend_list = legends,
    merge_legend = TRUE,
    heatmap_legend_side    = heatmap_legend_side,
    annotation_legend_side = annotation_legend_side,
    padding = unit(c(2, 2, 2, 2), "mm")
  )
  dev.off()
}

blend_white <- function(col, p = 0.55) {
  scales::colour_ramp(c("white", col))(p)
}

class_cols <- c(
  "PC"     = "#9FB7B5",  
  "SM"     = "#B9B7CC",  
  "LPC"    = "#A9BED6",  
  "PC O-"  = "#D9B0B8", 
  "TG"     = "#E6D8A8",
  "Others" = "#E2E2E2"   
)

# ------------------------------------------------------------------------------
# 3) LOAD INPUTS
# ------------------------------------------------------------------------------

bundle <- readRDS(clust_path)

mat_plot  <- bundle$mat_plot   # expected: lipids x timepoints (rows = lipids)
row_split <- bundle$row_split  # cluster label per lipid (aligned to rows of mat_plot)
ann_vec   <- bundle$ann_vec    # Up/Down/NS per lipid (aligned to rows of mat_plot)
class_vec <- bundle$class_vec
K_FINAL   <- bundle$K_FINAL
tp_levels <- bundle$tp_levels

# Use CLIP_Z from the clustering bundle when present (same clipping as
# 04_lip_clustering_time_trajectories.R)
CLIP_Z_BUNDLE <- bundle$CLIP_Z
if (!is.null(CLIP_Z_BUNDLE) && is.finite(CLIP_Z_BUNDLE)) {
  CLIP_Z <- CLIP_Z_BUNDLE
  pal_z <- circlize::colorRamp2(
    c(-CLIP_Z, 0, CLIP_Z),
    c("#3B6FB6", "#F7F7F7", "#B33A3A")
  )
}

mat_plot <- as.matrix(mat_plot)
stopifnot(nrow(mat_plot) == length(row_split), nrow(mat_plot) == length(ann_vec))

# Ensure timepoint columns order
stopifnot(all(tp_levels %in% colnames(mat_plot)))
mat_plot <- mat_plot[, tp_levels, drop = FALSE]

# Standardize labels
row_split <- factor(as.character(row_split), levels = as.character(seq_len(K_FINAL)))

ann_top <- factor(
  dplyr::recode(as.character(ann_vec),
                "Down" = "Down at TEA",
                "Up"   = "Up at TEA",
                "NS"   = "Not significant"),
  levels = c("Down at TEA", "Not significant", "Up at TEA")
)

stats <- read_delim(stats_path, delim = "\t", show_col_types = FALSE)

# ------------------------------------------------------------------------------
# PANEL A – Temporal heatmap
# ------------------------------------------------------------------------------

# A1) Transpose: rows=timepoints, cols=lipids
mat_t <- t(mat_plot)  
stopifnot(identical(rownames(mat_t), tp_levels))

# Column split (lipids) by cluster labels from bundle (aligned to columns of mat_t)
col_split <- row_split

# A2) Top annotation: T4 vs T0 direction (per lipid/column)
ha_top <- HeatmapAnnotation(
  `T4 vs T0` = ann_top,
  col = list(`T4 vs T0` = t4_cols),
  simple_anno_size = unit(ANNO_THICK_MM, "mm"),
  show_annotation_name = FALSE,
  which = "column"
)

# A3) Heatmap object
ht_A <- Heatmap(
  mat_t,
  name = "Z-score",
  col  = pal_z,
  
  cluster_rows = FALSE,
  row_order    = tp_levels,
  
  cluster_columns       = TRUE,
  column_order          = colnames(mat_t),
  column_split          = col_split,
  column_gap            = unit(COLUMN_GAP_MM, "mm"),
  column_title_gp       = gpar(fontsize = BASE_SIZE - 2),
  
  show_column_names = FALSE,
  show_row_names    = TRUE,
  row_names_gp      = gpar(fontsize = BASE_SIZE - 1),
  
  top_annotation = ha_top,
  border         = HEATMAP_BORDER,
  rect_gp        = gpar(col = NA),
  
  heatmap_legend_param = list(
    title = "Z-score",
    at = c(-CLIP_Z, 0, CLIP_Z),
    direction = "horizontal"
  )
)

save_ht_svg(
  ht_A, out_A,
  width = 4.8, height = 3.2,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom"
)

# ==============================================================================
# PANEL B – Cluster profile plots + class composition pies
# ==============================================================================

# Colours refer to the k-means cluster numbers from
# 04_lip_clustering_time_trajectories.R (fixed seed)
cluster_cols <- c(
  "2" = COL_DOWN,     
  "1" = "grey70",     
  "3" = "#C8A0A8",    
  "4" = COL_UP        
)

CLUSTER_ORDER_B <- c("2", "1", "3", "4")

# 1) Long table: lipid trajectories
mat_z <- as.matrix(mat_plot)[, tp_levels, drop = FALSE]
stopifnot(nrow(mat_z) == length(row_split))

trend_df <- mat_z %>%
  as.data.frame() %>%
  tibble::rownames_to_column("LipidName") %>%
  mutate(Cluster = factor(as.character(row_split), levels = CLUSTER_ORDER_B)) %>%
  tidyr::pivot_longer(all_of(tp_levels), names_to = "Timepoint", values_to = "z") %>%
  mutate(Timepoint = factor(Timepoint, levels = tp_levels)) %>%
  filter(!is.na(Cluster))

trend_sum <- trend_df %>%
  dplyr::group_by(Cluster, Timepoint) %>%
  dplyr::summarise(
    med = median(z, na.rm = TRUE),
    q1  = quantile(z, 0.25, na.rm = TRUE),
    q3  = quantile(z, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

lab_df <- dplyr::distinct(trend_sum, Cluster) %>%
  mutate(Timepoint = tp_levels[1], y = 1.65)

# 2) Profile plot (left)
p_prof <- ggplot() +
  geom_line(
    data = trend_df,
    aes(Timepoint, z, group = LipidName, color = Cluster),
    alpha = 0.18, linewidth = 0.32
  ) +
  geom_ribbon(
    data = trend_sum,
    aes(Timepoint, ymin = q1, ymax = q3, fill = Cluster),
    alpha = 0.18
  ) +
  geom_line(
    data = trend_sum,
    aes(Timepoint, med, group = Cluster, color = Cluster),
    linewidth = 2.2
  ) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.5, alpha = 0.5) +
  geom_text(
    data = lab_df,
    aes(Timepoint, y, label = as.character(Cluster), color = Cluster),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 5,
    hjust = 0
  ) +
  facet_wrap(~Cluster, ncol = 1) +
  coord_cartesian(ylim = c(-2, 2)) +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) +
  scale_color_manual(values = cluster_cols, drop = FALSE) +
  scale_fill_manual(values = cluster_cols, drop = FALSE) +
  labs(x = NULL, y = "Lipid abundance (z-score)") +
  theme_qomics() +
  theme(
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.background = element_blank(),
    strip.text       = element_blank(),
    legend.position  = "none"
  )

# 3) Pie charts (right): class composition per cluster
pie_df <- tibble::tibble(
  Cluster = factor(as.character(row_split), levels = CLUSTER_ORDER_B),
  Class   = as.character(class_vec)
) %>%
  mutate(
    Class = dplyr::if_else(is.na(Class) | Class == "" | !Class %in% names(class_cols),
                           "Others", Class)
  ) %>%
  count(Cluster, Class, name = "n") %>%
  group_by(Cluster) %>%
  mutate(
    prop      = as.numeric(n / sum(n)),                 
    perc_lab  = paste0(round(prop * 100, 1), "%"),
    n_cluster = sum(n)
  ) %>%
  ungroup() %>%
  arrange(Cluster, desc(prop))

p_pie <- ggplot(pie_df, aes(x = 1, y = prop, fill = Class)) +
  geom_col(width = 1, color = "white", linewidth = 0.6) +
  geom_text(
    aes(x = 1, label = perc_lab),                       
    position = position_stack(vjust = 0.5),
    color = "black", size = 3.4, fontface = "bold"
  ) +
  geom_text(
    aes(x = 1.6, label = Class),                        
    position = position_stack(vjust = 0.5),
    hjust = 0, size = 3.2
  ) +
  coord_polar(theta = "y", clip = "off") +
  scale_x_continuous(limits = c(0.5, 1.6)) +
  facet_wrap(
    ~Cluster, ncol = 2,
    labeller = labeller(Cluster = function(x) {
      n_lab <- pie_df$n_cluster[match(x, pie_df$Cluster)]
      paste0("Cluster ", x, "\n(n = ", n_lab, ")")
    })
  ) +
  scale_fill_manual(values = class_cols, drop = FALSE) +
  theme_void(base_size = BASE_SIZE) +
  theme(
    strip.text = element_text(size = BASE_SIZE - 2, hjust = 0.5),
    plot.margin = margin(6, 18, 6, 6),
    legend.position = "none"
  )

# 4) Combine + export
p_B <- p_prof + p_pie + patchwork::plot_layout(widths = c(1, 2.5))
print(p_B)

save_svg(p_B, out_B, width = 6.2, height = 5.0)

# ==============================================================================
# PANEL C – Volcano plot (T4 vs T0)
# ==============================================================================

N_LABEL_UP   <- 6
N_LABEL_DOWN <- 6

volcano_df <- stats %>%
  transmute(
    LipidName = str_trim(LipidName),
    beta      = beta_T4_vs_T0,
    FDR       = pmax(FDR_T4_T0, 1e-300)
  ) %>%
  mutate(
    neglog10_fdr = -log10(FDR),
    endpoint_class = case_when(
      FDR < FDR_CUTOFF & beta >  FC_CUTOFF ~ "Up",
      FDR < FDR_CUTOFF & beta < -FC_CUTOFF ~ "Down",
      TRUE ~ "NS"
    ),
    endpoint_class = factor(endpoint_class, levels = c("Down", "NS", "Up"))
  )

# Label selection: top N by FDR (tie-break |β|) in each direction 
pick_top_labels <- function(df, direction = c("Up", "Down"), n = 6) {
  direction <- match.arg(direction)
  df %>%
    filter(endpoint_class == direction) %>%
    arrange(FDR, desc(abs(beta))) %>%
    slice_head(n = min(n, nrow(.))) %>%
    pull(LipidName)
}

LABEL_UP   <- pick_top_labels(volcano_df, "Up",   n = N_LABEL_UP)
LABEL_DOWN <- pick_top_labels(volcano_df, "Down", n = N_LABEL_DOWN)
LABEL_LIPIDS <- unique(c(LABEL_UP, LABEL_DOWN))

volcano_df <- volcano_df %>%
  mutate(label = if_else(LipidName %in% LABEL_LIPIDS, LipidName, NA_character_))

message("[Volcano labels] Up:   ", paste(LABEL_UP, collapse = ", "))
message("[Volcano labels] Down: ", paste(LABEL_DOWN, collapse = ", "))

# Aesthetics driven by significance (cap at 99.5th percentile) 
cap_val <- as.numeric(quantile(volcano_df$neglog10_fdr, 0.995, na.rm = TRUE))

volcano_df <- volcano_df %>%
  mutate(
    sig_cap  = pmin(neglog10_fdr, cap_val),
    pt_size  = rescale(sig_cap, to = c(1.8, 3.4), from = range(sig_cap, na.rm = TRUE)),
    pt_alpha = rescale(sig_cap, to = c(0.25, 1), from = range(sig_cap, na.rm = TRUE))
  )

# Counts and thresholds
up_count   <- sum(volcano_df$endpoint_class == "Up",   na.rm = TRUE)
down_count <- sum(volcano_df$endpoint_class == "Down", na.rm = TRUE)

y_thr <- -log10(FDR_CUTOFF)
min_x <- min(volcano_df$beta, na.rm = TRUE)
max_x <- max(volcano_df$beta, na.rm = TRUE)

# Plot 
p_C <- ggplot(volcano_df, aes(x = beta, y = neglog10_fdr)) +
  # thresholds
  geom_hline(
    yintercept = y_thr,
    linetype   = "dashed",
    color      = alpha("grey55", 0.65),
    linewidth  = 0.4
  ) +
  geom_vline(
    xintercept = c(-FC_CUTOFF, FC_CUTOFF),
    linetype   = "dashed",
    color      = alpha("grey55", 0.65),
    linewidth  = 0.4
  ) +
  # points: NS first (background)
  geom_point(
    data = volcano_df %>% filter(endpoint_class == "NS"),
    aes(size = pt_size, alpha = pt_alpha),
    color = COL_NS
  ) +
  # points: significant (colored)
  geom_point(
    data = volcano_df %>% filter(endpoint_class != "NS"),
    aes(color = endpoint_class, size = pt_size, alpha = pt_alpha)
  ) +
  scale_color_manual(values = c("Up" = COL_UP, "Down" = COL_DOWN)) +
  scale_size_identity() +
  scale_alpha_identity() +
  # labels
  geom_text_repel(
    data = volcano_df %>% filter(!is.na(label)),
    aes(label = label),
    box.padding       = 0.55,
    point.padding     = 0.25,
    size              = 3.3,
    max.overlaps      = Inf,
    segment.color     = "grey70",
    min.segment.length = 0
  ) +
  annotate(
    "text", x = min_x - 0.9, y = y_thr + 0.2,
    label = paste0("▼ ", down_count),
    color = COL_DOWN, fontface = "bold", size = 5.0, hjust = 0
  ) +
  annotate(
    "text", x = max_x + 0.25, y = y_thr + 0.2,
    label = paste0("▲ ", up_count),
    color = COL_UP, fontface = "bold", size = 5.0, hjust = 1
  ) +
  labs(
    x = expression("Effect size (T4 vs T0)"),
    y = expression(-log[10] * "(FDR)")
  ) +
  coord_cartesian(xlim = c(min_x - 0.7, max_x + 0.1)) +
  theme_qomics() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

print(p_C)
save_svg(p_C, out_C, width = 4.0, height = 4.0)

# ==============================================================================
# PANEL D – LION enrichment barplot 
# ==============================================================================

lion_raw <- read_delim(lion_path, delim = ",", show_col_types = FALSE)

MIN_ANNOTATED  <- 10
N_PER_DIR      <- 10
MAX_PER_FAMILY <- 2

# helper: colorbar guide (same settings)
cb_guide <- function() {
  guide_colorbar(
    title.position = "left",
    title.theme    = element_text(angle = 90, vjust = 0.5, size = BASE_SIZE - 1),
    barwidth       = unit(0.35, "cm"),
    barheight      = unit(2.4, "cm")
  )
}

# clean + filter 
lion_f <- lion_raw %>%
  mutate(
    FDR       = as.numeric(FDR),
    ES        = as.numeric(ES),
    Regulated = as.character(Regulated),
    
    Description = str_replace(Description, "^([a-z])", toupper),
    Description = str_replace(
      Description,
      "Ceramide phosphocholines \\(sphingomyelins\\)",
      "Sphingomyelins"
    ),
    
    Family = case_when(
      str_detect(tolower(Description), "triacyl|triradyl|lipid storage|lipid droplet") ~ "Storage/TG",
      str_detect(tolower(Description), "ether|alkyl,2-acyl") ~ "Ether lipids",
      str_detect(tolower(Description), "sphing|sphingomyel") ~ "Sphingolipids/SM",
      str_detect(tolower(Description), "plasma membrane|membrane component|endosome/lysosome") ~ "Compartment",
      str_detect(tolower(Description), "bilayer|lateral diffusion|transition temperature|intrinsic curvature") ~ "Biophysics",
      str_detect(tolower(Description), "headgroup") ~ "Headgroup",
      str_detect(tolower(Description), "glycerophospho|diacylglycerophospho|monoacylglycerophospho|\\[gp") ~ "Phospholipids (GP/PC)",
      str_detect(tolower(Description), "fatty acid|c\\d{1,2}:\\d") ~ "Fatty-acid features",
      TRUE ~ "Other"
    ),
    
    neglog10FDR = -log10(pmax(FDR, 1e-300))
  ) %>%
  filter(
    !is.na(FDR), !is.na(ES), !is.na(Regulated),
    FDR < FDR_CUTOFF,
    Annotated >= MIN_ANNOTATED
  )

# picker
pick_top <- function(df, n, max_per_family) {
  df_ord <- df %>% arrange(FDR, desc(abs(ES)))
  
  picked <- df_ord %>%
    group_by(Family) %>%
    slice_head(n = max_per_family) %>%
    ungroup() %>%
    arrange(FDR, desc(abs(ES))) %>%
    slice_head(n = min(n, nrow(.)))
  
  if (nrow(picked) < n) {
    picked <- bind_rows(
      picked,
      df_ord %>%
        anti_join(picked, by = "ID") %>%
        slice_head(n = n - nrow(picked))
    )
  }
  
  picked %>% slice_head(n = min(n, nrow(.)))
}

lion_sel <- bind_rows(
  pick_top(filter(lion_f, Regulated == "UP"),   n = N_PER_DIR, max_per_family = MAX_PER_FAMILY),
  pick_top(filter(lion_f, Regulated == "DOWN"), n = N_PER_DIR, max_per_family = MAX_PER_FAMILY)
) %>%
  mutate(Direction = factor(Regulated, levels = c("UP", "DOWN")))

# ordering
ord_desc <- lion_sel %>%
  arrange(Direction, desc(ES)) %>%
  pull(Description) %>%
  unique()

lion_sel <- lion_sel %>%
  mutate(Description = factor(Description, levels = ord_desc))

cap_fdr <- quantile(lion_sel$neglog10FDR, 0.99, na.rm = TRUE)
lion_sel <- lion_sel %>% mutate(neglog10_cap = pmin(neglog10FDR, cap_fdr))

# plot
p_D <- ggplot(lion_sel, aes(x = ES, y = Description)) +
  geom_vline(xintercept = 0, color = "grey30", linewidth = 0.5) +
  
  geom_segment(
    aes(x = 0, xend = ES, yend = Description),
    color = "grey60",
    linewidth = 1.0
  ) +
  
  # UP points + UP fill scale
  geom_point(
    data = filter(lion_sel, Regulated == "UP"),
    aes(size = Annotated, fill = neglog10_cap),
    shape = 21, color = "black", stroke = 0.25
  ) +
  scale_fill_gradient(
    name  = expression(-log[10] * "(FDR)"),
    low   = "#F6DADA",
    high  = COL_UP,
    guide = cb_guide()
  ) +
  
  ggnewscale::new_scale_fill() +
  
  # DOWN points + DOWN fill scale
  geom_point(
    data = filter(lion_sel, Regulated == "DOWN"),
    aes(size = Annotated, fill = neglog10_cap),
    shape = 21, color = "black", stroke = 0.25
  ) +
  scale_fill_gradient(
    name  = expression(-log[10] * "(FDR)"),
    low   = "#DCEEEE",
    high  = COL_DOWN,
    guide = cb_guide()
  ) +
  
  scale_size_continuous(
    name   = "Lipid count",
    range  = c(3, 7),
    breaks = pretty_breaks(n = 4),
    guide  = guide_legend(
      title.position = "left",
      title.theme = element_text(angle = 90, vjust = 0.5, size = BASE_SIZE - 1),
      keyheight = unit(0.35, "cm"),
      keywidth  = unit(0.35, "cm")
    )
  ) +
  
  scale_y_discrete(limits = rev(levels(lion_sel$Description))) +
  labs(x = "Enrichment score", y = NULL) +
  theme_qomics() +
  theme(
    axis.text.y        = element_text(size = BASE_SIZE - 1),
    axis.title.x       = element_text(margin = margin(t = 10)),
    legend.position    = "right",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank()
  )

print(p_D)
save_svg(p_D, out_D, width = 7.2, height = 5.2)

# ==============================================================================
# PANEL E – Acyl-chain structural space for membrane lipid classes
# ==============================================================================

KEEP_CLASSES <- c("PC", "PC O-", "SM")
PAD_X <- 1
PAD_Y <- 1

# Parse lipid shorthand -> Class + CN + DB
get_class <- function(x) {
  x <- str_trim(x)
  case_when(
    str_detect(x, "^PC\\s+O-") ~ "PC O-",
    str_detect(x, "^PC\\s+P-") ~ "PC P-",
    str_detect(x, "^PC\\s+")   ~ "PC",
    str_detect(x, "^LPC\\s+")  ~ "LPC",
    str_detect(x, "^SM\\s+")   ~ "SM",
    str_detect(x, "^TG\\s+")   ~ "TG",
    TRUE ~ NA_character_
  )
}

# Extract first CN:DB occurrence (works for "PC 36:5" and "SM 28:1;O2")
get_cn_db_total <- function(x) {
  comp <- str_extract(x, "\\d+[:]\\d+")
  if (is.na(comp)) return(tibble(CN = NA_integer_, DB = NA_integer_))
  parts <- str_split(comp, ":", simplify = TRUE)
  tibble(CN = as.integer(parts[, 1]), DB = as.integer(parts[, 2]))
}

# Sum chains for TG ("TG 16:0_18:2_18:3")
get_cn_db_sum <- function(x) {
  rhs <- str_remove(x, "^TG\\s+")
  chains <- unlist(str_split(rhs, "[_/]"))
  chains <- chains[str_detect(chains, "^\\d+[:]\\d+$")]
  if (length(chains) == 0) return(tibble(CN = NA_integer_, DB = NA_integer_))
  cn <- sum(as.integer(str_extract(chains, "^\\d+")))
  db <- sum(as.integer(str_extract(chains, "(?<=:)\\d+")))
  tibble(CN = cn, DB = db)
}

parse_lipid_one <- function(x) {
  cls <- get_class(x)
  if (is.na(cls)) return(tibble(Class = NA_character_, CN = NA_integer_, DB = NA_integer_))
  cndb <- if (cls == "TG") get_cn_db_sum(x) else get_cn_db_total(x)
  tibble(Class = cls, CN = cndb$CN, DB = cndb$DB)
}

parse_lipid <- function(vec) bind_rows(lapply(vec, parse_lipid_one))

plot_df <- stats %>%
  transmute(
    LipidName = LipidName,
    effect    = beta_T4_vs_T0,
    FDR       = FDR_T4_T0
  ) %>%
  bind_cols(parse_lipid(.$LipidName)) %>%
  filter(Class %in% KEEP_CLASSES) %>%
  filter(!is.na(CN), !is.na(DB)) %>%
  mutate(
    sig = !is.na(FDR) & FDR < FDR_CUTOFF,
    neglog10FDR = -log10(pmax(FDR, 1e-300)),
    Class = factor(Class, levels = KEEP_CLASSES)
  )

# Zoom window per class driven by significant lipids (plus padding)
rng_tbl <- plot_df %>%
  filter(sig) %>%
  group_by(Class) %>%
  summarise(
    xmin = min(CN, na.rm = TRUE) - PAD_X,
    xmax = max(CN, na.rm = TRUE) + PAD_X,
    ymin = min(DB, na.rm = TRUE) - PAD_Y,
    ymax = max(DB, na.rm = TRUE) + PAD_Y,
    .groups = "drop"
  )

plot_df <- plot_df %>% filter(Class %in% rng_tbl$Class)

plot_zoom <- plot_df %>%
  left_join(rng_tbl, by = "Class") %>%
  filter(CN >= xmin, CN <= xmax, DB >= ymin, DB <= ymax)

# Symmetric limits for effect size (balanced diverging palette)
eff_lim <- max(abs(plot_zoom$effect), na.rm = TRUE)

p_E <- ggplot(plot_zoom, aes(CN, DB)) +
  geom_point(
    data = plot_zoom %>% filter(!sig),
    color = "grey70", size = 1.2, alpha = 0.55
  ) +
  geom_point(
    data = plot_zoom %>% filter(sig),
    aes(size = neglog10FDR, fill = effect),
    shape = 21, color = "black", stroke = 0.5
  ) +
  facet_wrap(~ Class, scales = "free") +
  scale_fill_gradient2(
    name = "Effect size\n(T4 vs T0)",
    low = COL_DOWN,
    mid = "white",
    high = COL_UP,
    midpoint = 0,
    limits = c(-eff_lim, eff_lim),
    oob = squish,
    breaks = pretty_breaks(n = 5),
    guide = guide_colourbar(
      title.position = "left",
      title.theme = element_text(
        angle = 90,
        vjust = 0.5,
        size  = BASE_SIZE - 2
      ),
      barwidth  = unit(0.45, "cm"),
      barheight = unit(2.4, "cm")
    )
  ) +
  scale_size_continuous(
    name = expression(-log[10] * "(FDR)"),
    range = c(2.3, 7.5),
    breaks = pretty_breaks(n = 4),
    guide = guide_legend(
      title.position = "left",
      title.theme = element_text(
        angle = 90,
        vjust = 0.5,
        size  = BASE_SIZE - 2
      ),
      keyheight = unit(0.45, "cm"),
      keywidth  = unit(0.45, "cm")
    )
  ) +
  scale_x_continuous(breaks = breaks_width(2), minor_breaks = breaks_width(1)) +
  scale_y_continuous(breaks = breaks_width(1), minor_breaks = breaks_width(0.5)) +
  labs(
    x = "Number of carbon atoms in side chain",
    y = "Saturation"
  ) +
  theme_qomics() +
  theme(
    strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.3),
    strip.text = element_text(face = "bold", size = BASE_SIZE),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.35),
    panel.grid.minor = element_line(color = "grey93", linewidth = 0.25),
    legend.key.height = unit(0.55, "cm"),
    legend.key.width  = unit(0.45, "cm")
  )

print(p_E)
save_svg(p_E, out_E, width = 11.5, height = 3.5)

# ------------------------------------------------------------------------------
# DONE
# ------------------------------------------------------------------------------
cat(
  "\n[OK] Figure 5 panels saved:\n",
  " - ", out_A, "\n",
  " - ", out_B, "\n",
  " - ", out_C, "\n",
  " - ", out_D, "\n",
  " - ", out_E, "\n",
  sep = ""
)
