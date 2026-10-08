# ==============================================================================
# Project: Preterm EV Multi-Omics
# Figure:  Figure 3 – Datasets overview
# Script:  fig3_dataset_overview.R
# Author:  Nicole Grinovero
# Date:    Jan 2026
#
# Panels:
#   a) Proteomics dynamic range (rank plot)
#   b) Proteomics subcellular localization (Enrichr TOP10)
#   c) Lipidomics class composition
#
# Inputs:
#   - data/proteomics/processed/proteins_norm_imputed.txt
#   - outputs/tables/preterm/proteomics/enrich_Jensen_COMPARTMENTS.txt
#   - data/lipidomics/processed/lipids_cleaned.txt
#
# Outputs:
#   - outputs/figures/fig3_dataset_overview/fig3a_proteomics_rankplot.svg
#   - outputs/figures/fig3_dataset_overview/fig3b_proteomics_compartments.svg
#   - outputs/figures/fig3_dataset_overview/fig3c_lipidomics_class_barplot.svg
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
  library(cowplot)
  library(svglite)
})

# ------------------------------------------------------------------------------
# 0) Global style
# ------------------------------------------------------------------------------
BASE_SIZE <- 12

COL_ACCENT <- "#7F4E57"
COL_ACCENT_FILL <- "#B88A92"
COL_BAR <- "#9FAAA9"
COL_HL <- "#6C8C8A"
COL_EDGE <- "#6E6E6E"

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
out_dir <- here("outputs", "figures", "fig3_dataset_overview")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

prot_file <- here("data", "proteomics", "processed", "proteins_norm_imputed.txt")
lip_file <- here("data", "lipidomics", "processed", "lipids_cleaned.txt")
enrichr_top10_file <- here(
  "outputs", "tables", "preterm", "proteomics",
  "enrich_Jensen_COMPARTMENTS.txt"
)

# ------------------------------------------------------------------------------
# PANEL A – Proteomics rank plot (dynamic range)
# ------------------------------------------------------------------------------
marker_genes <- c(
  "APOA1", "ALB", "LGALS3BP", "NCAM1", "SPP1",
  "SDCBP", "STX7", "STEAP3", "ADAM10", "SNX3",
  "APP", "L1CAM"
)

expr <- read_delim(prot_file, delim = "\t", show_col_types = FALSE) %>%
  mutate(GeneName = na_if(str_squish(GeneName), ""))

quant_cols <- names(expr)[vapply(expr, is.numeric, logical(1))]
stopifnot(length(quant_cols) > 0)

rank_df <- expr %>%
  mutate(Median = apply(select(., all_of(quant_cols)), 1, median, na.rm = TRUE)) %>%
  filter(!is.na(GeneName)) %>%
  group_by(GeneName) %>%
  slice_max(Median, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(Median)) %>%
  mutate(
    Rank   = row_number(),
    Marker = if_else(GeneName %in% marker_genes, "Marker", "Other")
  )

p_b <- ggplot(rank_df, aes(x = Rank, y = Median)) +
  geom_point(color = "grey55", alpha = 0.22, size = 1.7) +
  geom_point(
    data = filter(rank_df, Marker == "Marker"),
    shape = 21,
    fill = COL_ACCENT_FILL,
    color = COL_ACCENT,
    stroke = 0.7,
    size = 2.7,
    alpha = 0.95
  ) +
  ggrepel::geom_text_repel(
    data = filter(rank_df, Marker == "Marker"),
    aes(label = GeneName),
    size = 3.5,
    color = "grey10",
    segment.size = 0.35,
    segment.color = COL_ACCENT,
    box.padding = 0.35,
    point.padding = 0.20,
    max.overlaps = Inf,
    seed = 123
  ) +
  labs(
    x = "Protein rank",
    y = expression("Median log"[2] * "-abundance")
  ) +
  theme_qomics() +
  theme(panel.grid = element_blank())

print(p_b)

save_svg(
  plot     = p_b,
  filename = file.path(out_dir, "fig3a_proteomics_rankplot.svg"),
  width    = 3.4,
  height   = 2.7
)

# ------------------------------------------------------------------------------
# PANEL B – Proteomics subcellular localization (Enrichr TOP10 table)
# ------------------------------------------------------------------------------
df_enr <- read_delim(enrichr_top10_file, delim = "\t", show_col_types = FALSE)

req_enr <- c("Term", "Count", "Odds.Ratio", "Adjusted.P.value")
stopifnot(all(req_enr %in% names(df_enr)))

df_enr <- df_enr %>%
  mutate(
    neglog10_adjP = -log10(Adjusted.P.value),
    Term = str_replace_all(Term, "_", " ")
  ) %>%
  arrange(desc(Count)) %>%
  slice_head(n = 10) %>%
  mutate(Term = factor(Term, levels = rev(unique(Term))))

p_c <- ggplot(df_enr, aes(x = Count, y = Term)) +
  geom_point(aes(size = neglog10_adjP, color = Odds.Ratio), alpha = 0.90) +
  scale_color_gradient(
    low  = "#C8A0A8",
    high = "#7F4E57",
    name = "Odds ratio"
  ) +
  scale_size(
    range = c(2, 6.5),
    name  = expression(-log[10] ~ "(adj. P)")
  ) +
  guides(
    color = guide_colorbar(
      title.position = "left",
      title.theme = element_text(angle = 90, vjust = 0.5, size = BASE_SIZE - 1),
      barwidth = unit(0.35, "cm"),
      barheight = unit(2.2, "cm")
    ),
    size = guide_legend(
      title.position = "left",
      title.theme = element_text(angle = 90, vjust = 0.5)
    )
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(x = "Protein count", y = NULL) +
  theme_qomics() +
  theme(
    axis.text.y = element_text(size = 10.5),
    panel.grid.major.x = element_line(linetype = "dashed", linewidth = 0.3, color = "grey80"),
    panel.grid.major.y = element_line(linewidth = 0.25, color = "grey92"),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = BASE_SIZE - 2),
    legend.text = element_text(size = BASE_SIZE - 3),
    legend.key.height = unit(0.35, "cm"),
    legend.key.width = unit(0.35, "cm"),
    legend.spacing.y = unit(0.05, "cm"),
    legend.box.spacing = unit(0.15, "cm")
  )

print(p_c)

save_svg(
  plot     = p_c,
  filename = file.path(out_dir, "fig3b_proteomics_compartments.svg"),
  width    = 4.6,
  height   = 2.8
)

# ------------------------------------------------------------------------------
# PANEL C – Lipidomics class composition
# ------------------------------------------------------------------------------
highlight_classes <- c("SM")
MIN_CLASS_COUNT <- 2
RARE_LABEL <- "Other"

facet_labels <- c(
  "Glycerophospholipids" = "Glycerophospholipids",
  "Glycerolipids"        = "Glycerolipids",
  "Sphingolipids"        = "Sphingolipids",
  "Fatty acyls"          = "Fatty\nacyls",
  "Sterol lipids"        = "Sterol\nlipids"
)

category_order <- names(facet_labels)

lip <- read_delim(lip_file, delim = "\t", show_col_types = FALSE) %>% 
  filter(!str_detect(LipidName, "\\(d[0-9]+\\)"))

stopifnot(all(c("Category", "Class") %in% names(lip)))

class_df <- lip %>%
  count(Category, Class, name = "n") %>%
  mutate(Class = if_else(n < MIN_CLASS_COUNT, RARE_LABEL, Class)) %>%
  group_by(Category, Class) %>%
  summarise(n = sum(n), .groups = "drop")

plot_df <- class_df %>%
  mutate(
    Category  = factor(Category, levels = category_order),
    Highlight = if_else(Class %in% highlight_classes, "Highlight", "Other"),
    Highlight = if_else(Class == RARE_LABEL, "Other", Highlight)
  ) %>%
  group_by(Category) %>%
  mutate(
    Class = fct_reorder(Class, n, .desc = TRUE),
    Class = fct_relevel(Class, RARE_LABEL, after = Inf)
  ) %>%
  ungroup()

p_d <- ggplot(plot_df, aes(x = Class, y = n)) +
  geom_col(
    aes(fill = Highlight),
    width = 0.85,
    color = COL_EDGE,
    linewidth = 0.18
  ) +
  geom_text(
    aes(label = n),
    vjust = -0.25,
    size = 3.5,
    color = "black"
  ) +
  facet_grid(
    ~Category,
    scales = "free_x",
    space = "free_x",
    labeller = labeller(Category = facet_labels)
  ) +
  scale_fill_manual(values = c("Other" = COL_BAR, "Highlight" = COL_HL)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.14))) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Number of lipid species") +
  theme_qomics() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = BASE_SIZE - 1, lineheight = 0.95),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

print(p_d)

save_svg(
  plot     = p_d,
  filename = file.path(out_dir, "fig3c_lipidomics_class_barplot.svg"),
  width    = 8,
  height   = 3
)

cat("\nFigure 3 panels saved in:\n", out_dir, "\n", sep = "")
