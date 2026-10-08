# ==============================================================================
# Fig. S1 — Proteomics technical validation
# Project:   Plasma EV Multi-Omics — Technical validation
# Author:    Nicole Grinovero
# Date:      Feb 2026
#
# Output:
#   - Fig. S1a: Sample-to-sample correlation heatmap (log10 intensity; Pearson R)
#   - Fig. S1b: Protein group overlap across conditions (UpSet plot)
#
# Inputs:
#   - data/techrep/proteomics/20260206_182115_plasma_EV_multiomics_techrep_Report.tsv
#     (Spectronaut report; *.PG.Quantity columns)
#
# Outputs:
#   outputs/figures/figS1/
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(stringr)
  library(ComplexHeatmap)
  library(circlize)
  library(ComplexUpset)
  library(svglite)
  library(scales)
  library(grid)
  library(ragg)
})

# ------------------------------------------------------------------------------
# 0) PATHS
# ------------------------------------------------------------------------------
report_path <- here("data","techrep","proteomics","20260206_182115_plasma_EV_multiomics_techrep_Report.tsv")

out_dir <- here("outputs","figures","figS1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1) SETTINGS
# ------------------------------------------------------------------------------
# Condition order (global, used in heatmap annotation + upset)
ord_levels <- c(
  "Mag-Net BAMM 10µL", "Mag-Net BAMM 5µL",
  "Mag-Net IPA 10µL",  "Mag-Net IPA 5µL",
  "Mag-Net MeOH 10µL", "Mag-Net MeOH 5µL",
  "Mag-Net standard 10µL","Mag-Net standard 5µL",
  "Neat"
)

# Condition palette (same as Fig. 1)
cond_colors <- c(
  "Neat"                  = "whitesmoke",
  "Mag-Net standard 5µL"  = "#B8B8B8",
  "Mag-Net standard 10µL" = "#9E9E9E",
  "Mag-Net MeOH 5µL"      = "#F3E8C9",
  "Mag-Net MeOH 10µL"     = "#E6D3A3",
  "Mag-Net BAMM 5µL"      = "#C7DDE8",
  "Mag-Net BAMM 10µL"     = "#9BBFD4",
  "Mag-Net IPA 5µL"       = "#E1D7E6",
  "Mag-Net IPA 10µL"      = "#C3B2D1"
)

# Robust µL normalization
norm_ul <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("μL", "µL") %>%
    str_replace_all("(?i)uL", "µL") %>%
    str_replace_all("\\s+", " ") %>%
    str_squish()
}

# Explicit recode map (Spectronaut condition tokens -> display labels)
recode_condition <- function(x) {
  out <- case_when(
    x == "Neat" ~ "Neat",
    x == "MagNet_standard_5uL"  ~ "Mag-Net standard 5µL",
    x == "MagNet_standard_10uL" ~ "Mag-Net standard 10µL",
    x == "MagNet_MeOH_5uL"      ~ "Mag-Net MeOH 5µL",
    x == "MagNet_MeOH_10uL"     ~ "Mag-Net MeOH 10µL",
    x == "MagNet_BAMM_5uL"      ~ "Mag-Net BAMM 5µL",
    x == "MagNet_BAMM_10uL"     ~ "Mag-Net BAMM 10µL",
    x == "MagNet_IPA_5uL"       ~ "Mag-Net IPA 5µL",
    x == "MagNet_IPA_10uL"      ~ "Mag-Net IPA 10µL",
    TRUE ~ as.character(x)
  )
  norm_ul(out)
}

is_det <- function(x) is.finite(x) & !is.na(x) & x > 0

# Save helpers
save_gg <- function(p, prefix, w, h, dpi = 320) {
  print(p)
  ggsave(paste0(prefix, ".svg"), p, width = w, height = h, device = svglite::svglite)
  ggsave(paste0(prefix, ".png"), p, width = w, height = h, dpi = dpi, bg = "white")
}

save_ht <- function(ht, prefix, w, h, res = 320) {
  # show in console
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  
  # SVG
  svglite::svglite(paste0(prefix, ".svg"), width = w, height = h)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  
  # PNG (high quality)
  ragg::agg_png(paste0(prefix, ".png"), width = w, height = h, units = "in", res = res, background = "white")
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
}

# ------------------------------------------------------------------------------
# 2) LOAD + PARSE REPORT.TSV
# ------------------------------------------------------------------------------
rep <- readr::read_delim(report_path, delim = "\t", show_col_types = FALSE)
stopifnot("PG.ProteinGroups" %in% colnames(rep))

sample_cols_raw <- grep("\\.raw\\.PG\\.Quantity$", colnames(rep), value = TRUE)
sample_cols_lin <- grep("\\.PG\\.Quantity$", colnames(rep), value = TRUE)
sample_cols <- if (length(sample_cols_raw) > 0) sample_cols_raw else sample_cols_lin

sample_meta <- tibble(sample_col = sample_cols) %>%
  mutate(
    sample_clean  = str_remove(sample_col, "^\\[\\d+\\]\\s+"),
    replicate     = str_extract(sample_clean, "C\\d+"),
    # strip the tail (works for both raw.PG.Quantity and PG.Quantity)
    condition_raw = sample_clean %>%
      str_remove("_C\\d+_\\d+\\.(raw\\.)?PG\\.Quantity$"),
    Condition     = recode_condition(condition_raw),
    Condition_ord = factor(Condition, levels = ord_levels)
  )

# long expression table
expr_long <- rep %>%
  select(PG.ProteinGroups, all_of(sample_cols)) %>%
  pivot_longer(all_of(sample_cols), names_to = "sample_col", values_to = "Intensity") %>%
  left_join(sample_meta, by = "sample_col") %>%
  mutate(Present = is_det(Intensity))

# ------------------------------------------------------------------------------
# 3) FIG S1A — CORRELATION HEATMAP (ComplexHeatmap)
# ------------------------------------------------------------------------------
expr_wide <- expr_long %>%
  select(PG.ProteinGroups, sample_col, Intensity) %>%
  pivot_wider(names_from = sample_col, values_from = Intensity)

mat <- expr_wide %>%
  select(-PG.ProteinGroups) %>%
  as.matrix()

mat_log <- log10(mat + 1)
cor_mat <- cor(mat_log, use = "pairwise.complete.obs", method = "pearson")

stopifnot(identical(colnames(cor_mat), sample_meta$sample_col))

ha <- HeatmapAnnotation(
  Condition = sample_meta$Condition_ord,
  col = list(Condition = cond_colors),
  gp = gpar(col = "black", lwd = 0.6),
  annotation_name_side = "left",
  annotation_legend_param = list(
    Condition = list(
      nrow = 9,
      border = TRUE,
      border_gp = gpar(col = "black", lwd = 0.6)
    )
  )
)

col_fun <- colorRamp2(
  c(0.85, 0.90, 0.95, 0.98, 1.00),
  c("#F2F2F2", "#D9E6F2", "#B5D3E8", "#6FAED6", "#2C7FB8")
)

ht <- Heatmap(
  cor_mat,
  name = "Pearson r",
  col = col_fun,
  border_gp = gpar(color = "black"),
  top_annotation = ha,
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_dend_width = unit(0.8, "cm"),
  column_dend_height = unit(0.8, "cm"),
  row_title = NULL,
  column_title = expression("Sample-to-sample correlation (" * log[10] * " intensity)"),
  heatmap_legend_param = list(
    title = "Correlation",
    at = c(0.85, 0.90, 0.95, 1.0),
    legend_gp = gpar(col = "black", lwd = 0.6)
  )
)

ht

# Save immediately
save_ht(
  ht,
  prefix = file.path(out_dir, "figS1a_proteomics_correlation_heatmap"),
  w = 6.5, h = 4.0
)

# ------------------------------------------------------------------------------
# 4) FIG S1B — UPSET OVERLAP (Protein groups across conditions)
# ------------------------------------------------------------------------------
# presence matrix: one row per PG, one logical col per condition
presence <- expr_long %>%
  filter(!is.na(Condition)) %>%
  group_by(Condition, PG.ProteinGroups) %>%
  summarise(Present = any(Present), .groups = "drop") %>%
  pivot_wider(names_from = Condition, values_from = Present, values_fill = FALSE)

# enforce global condition order and keep only existing columns
set_levels <- ord_levels[ord_levels %in% colnames(presence)]
presence <- presence %>% select(PG.ProteinGroups, all_of(set_levels))

# queries for coloring matrix dots by condition
q_list <- purrr::map(
  set_levels,
  ~ ComplexUpset::upset_query(set = .x, fill = cond_colors[[.x]])
)

p_upset <- ComplexUpset::upset(
  presence,
  intersect = set_levels,
  name = NULL,
  queries = q_list,
  min_size = 15,
  width_ratio = 0.25,
  sort_sets = "descending",
  base_annotations = list(
    "Intersection size" = (
      ComplexUpset::intersection_size(
        bar_number_threshold = 1,
        width = 0.75
      ) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black")
        )
    )
  ),
  stripes = ComplexUpset::upset_stripes(
    geom = geom_segment(linewidth = 7.5),
    colors = c("grey95", "white")
  ),
  matrix = ComplexUpset::intersection_matrix(
    geom = geom_point(
      shape = 21,
      size = 3.3,
      stroke = 0.3,
      colour = "black"
    ),
    geom_segment(linewidth = 0.8, colour = "grey35")
  ),
  set_sizes = (
    upset_set_size(geom = geom_bar(width = 0.70, color = "black", linewidth = 0.25)) +
      scale_fill_manual(values = cond_colors, breaks = set_levels, limits = set_levels, drop = FALSE) +
      geom_text(
        aes(label = after_stat(count)),
        stat = "count",
        vjust = 0.4, hjust = -0.25,
        size = 3.0,
        color = "black",
        fontface = "bold"
      ) +
      theme(
        panel.grid = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.ticks.x = element_line()
      )
  )
) +
  ggtitle("Protein group overlap across conditions") +
  theme(
    text = element_text(size = 13),
    plot.title = element_text(size = 12),
    panel.grid = element_blank(),
    axis.title.x = element_blank()
  )

# Save immediately (and show)
save_gg(
  p_upset,
  prefix = file.path(out_dir, "figS1b_proteomics_upset_overlap"),
  w = 7, h = 5.5
)

message("Done. Outputs written to: ", out_dir)

# Source data for Fig. S1b: presence of each protein group by condition
presence_annot <- presence %>%
  left_join(rep %>% select(PG.ProteinGroups, PG.Genes, PG.ProteinDescriptions), by = "PG.ProteinGroups")
write_delim(presence_annot, file.path(out_dir, "figS1b_protein_presence_by_condition.tsv"), delim = "\t")
