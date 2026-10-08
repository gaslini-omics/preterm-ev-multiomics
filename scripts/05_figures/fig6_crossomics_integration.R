# ==============================================================================
# Project: Preterm EV Multi-Omics – Integration (proteomics + lipidomics)
# Figure:  Figure 6 – Cross-omics integration 
# Script:  fig6_crossomics_integration.R
# Author:  Nicole Grinovero
# Date:    Jan 2026
#
# Panels:
#   B) Global landscape of protein–lipid correlations 
#   C) Cross-omics module composition (proteins vs lipids)
#   D) Functional interpretation of modules (GO:BP enrichment)
#
# Inputs:
#   - outputs/tables/preterm/integration/integration_correlation_landscape_stats.tsv
#   - outputs/tables/preterm/integration/integration_edges.tsv
#   - outputs/tables/preterm/integration/integration_modules_summary.tsv
#   - outputs/tables/preterm/integration/integration_modules_enrichment_GOBP.tsv
#
# Outputs:
#   - outputs/figures/fig6_crossomics_integration/fig6b_correlation_landscape.svg
#   - outputs/figures/fig6_crossomics_integration/fig6c_module_composition.svg
#   - outputs/figures/fig6_crossomics_integration/fig6d_module_enrichment_GOBP.svg
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(svglite)
})

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
# Palette / constants
# ------------------------------------------------------------------------------
COL_POS <- "#B33A3A"
COL_NEG <- "#3B6FB6"

# Composition colors
COL_PROT <- "#b9d5e1ff"
COL_LIP  <- "#f3cfadff"

MOD_LEVELS <- c("M1","M2","M3","M4","M5")

# module colors for Panel D 
MOD_COLS <- c(
  "M1" = "#C8D7C1",
  "M2" = "#F1D7C2",
  "M3" = "#C9BEDD",
  "M4" = "#D7D7D7",
  "M5" = "#BFD3DE"
)

# ------------------------------------------------------------------------------
# Helper
# ------------------------------------------------------------------------------
wrap_and_abbrev_go <- function(x, width = 18, max_lines = 2, max_chars = 60) {
  x2 <- stringr::str_squish(x)
  x2 <- stringr::str_replace(x2, "^([a-z])", toupper)
  
  wrapped <- stringr::str_wrap(x2, width = width)
  parts   <- strsplit(wrapped, "\n", fixed = TRUE)
  
  out <- vapply(parts, function(p) {
    # keep at most max_lines
    if (length(p) > max_lines) {
      p <- p[1:max_lines]
    }
    txt <- paste(p, collapse = "\n")
    
    # abbreviate ONLY if still too long
    if (nchar(gsub("\n", "", txt)) > max_chars) {
      txt <- paste0(substr(txt, 1, max_chars - 1), "…")
    }
    txt
  }, character(1))
  
  out
}

# ------------------------------------------------------------------------------
# 1) Paths
# ------------------------------------------------------------------------------
out_dir <- here("outputs", "figures", "fig6_crossomics_integration")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Integration inputs
stats_path <- here("outputs","tables","preterm","integration","integration_correlation_landscape_stats.tsv")
edges_path <- here("outputs","tables","preterm","integration","integration_edges.tsv")
sum_path   <- here("outputs","tables","preterm","integration","integration_modules_summary.tsv")
enr_path   <- here("outputs","tables","preterm","integration","integration_modules_enrichment_GOBP.tsv")

# Outputs
out_B <- file.path(out_dir, "fig6b_correlation_landscape.svg")
out_C <- file.path(out_dir, "fig6c_module_composition.svg")
out_D <- file.path(out_dir, "fig6d_module_enrichment_GOBP.svg")

# ==============================================================================
# PANEL B — Correlation landscape (ONLY edges used for the network)
# ==============================================================================
edges_df <- read_tsv(edges_path, show_col_types = FALSE) %>%
  mutate(
    rho     = as.numeric(rho),
    abs_rho = abs(rho),
    sign    = if_else(rho >= 0, "Positive", "Negative"),
    sign    = factor(sign, levels = c("Negative","Positive"))
  )

stats_tbl <- read_tsv(stats_path, show_col_types = FALSE)

n_edges <- nrow(edges_df)
n_pos   <- sum(edges_df$rho > 0, na.rm = TRUE)
n_neg   <- sum(edges_df$rho < 0, na.rm = TRUE)
med_abs <- median(edges_df$abs_rho, na.rm = TRUE)
rng_rho <- range(edges_df$rho, na.rm = TRUE)

ann_txt <- paste0(
  "Edges used for network\n",
  "FDR < 0.05 & |ρ| ≥ 0.35",
  "\n\nTotal tests: ", format(stats_tbl$n_tests_total, big.mark = ","),
  "\nEdges retained: ", format(n_edges, big.mark = ","),
  "\nNegative: ", format(n_neg, big.mark = ","),
  "\nPositive: ", format(n_pos, big.mark = ","),
  "\nMedian |ρ|: ", round(med_abs, 2),
  "\nρ range: [", round(rng_rho[1], 2), ", ", round(rng_rho[2], 2), "]"
)

p_B <- ggplot(edges_df, aes(x = rho, fill = sign)) +
  geom_histogram(
    bins = 60,
    color = "white",
    linewidth = 0.25,
    alpha = 0.60
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.7,
    color = "grey30"
  ) +
  scale_fill_manual(values = c("Positive" = COL_POS, "Negative" = COL_NEG)) +
  annotate(
    "label",
    x = -Inf, y = Inf,
    hjust = -0.02, vjust = 1.08,
    label = ann_txt,
    size = 2.8,
    label.padding = unit(0.35, "lines"),
    fill = alpha("white", 0.9),
    color = "grey15"
  ) +
  labs(
    x = "Spearman correlation coefficient (ρ)",
    y = "Number of protein–lipid pairs",
    fill = NULL
  ) +
  theme_qomics() +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.major.x = element_blank()
  )

print(p_B)
save_svg(p_B, out_B, width = 5.5, height = 4.2)

# ==============================================================================
# PANEL C — Module composition (Protein vs Lipid)
# ==============================================================================
# Module counts (integration_modules_summary.tsv, from 01_integrate_prot_lip_network.R)
mod_sum <- read_tsv(sum_path, show_col_types = FALSE) %>%
  transmute(ModuleLabel = as.character(ModuleLabel), n_Protein = Protein, n_Lipid = Lipid)
stopifnot(all(MOD_LEVELS %in% mod_sum$ModuleLabel))

mod_sum <- mod_sum %>%
  filter(ModuleLabel %in% MOD_LEVELS) %>%
  mutate(
    n_Total   = n_Protein + n_Lipid,
    frac_prot = n_Protein / n_Total,
    frac_lip  = n_Lipid   / n_Total,
    ModuleLabel = factor(ModuleLabel, levels = MOD_LEVELS)
  )

comp_long <- mod_sum %>%
  select(ModuleLabel, frac_prot, frac_lip) %>%
  pivot_longer(cols = starts_with("frac_"), names_to = "Layer", values_to = "Frac") %>%
  mutate(
    Layer = recode(Layer, frac_prot = "Protein", frac_lip = "Lipid"),
    Layer = factor(Layer, levels = c("Protein","Lipid"))
  )

y_top_n  <- 1.045
y_top_pl <- 1.000

p_C <- ggplot(comp_long, aes(x = ModuleLabel, y = Frac, fill = Layer)) +
  geom_col(width = 0.72, color = "grey25", linewidth = 0.25) +
  scale_fill_manual(
    values = c("Protein" = COL_PROT, "Lipid" = COL_LIP),
    name = NULL
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  geom_text(
    data = mod_sum,
    aes(x = ModuleLabel, y = y_top_n, label = paste0("n=", n_Total)),
    inherit.aes = FALSE,
    size = 3.2,
    color = "grey15"
  ) +
  geom_text(
    data = mod_sum,
    aes(x = ModuleLabel, y = y_top_pl, label = paste0("P=", n_Protein, "  L=", n_Lipid)),
    inherit.aes = FALSE,
    size = 3.0,
    color = "grey20"
  ) +
  coord_cartesian(ylim = c(0, 1.09), clip = "off") +
  labs(x = NULL, y = "Module composition") +
  theme_qomics() +
  theme(
    legend.position = "top",
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 12),
    plot.margin = margin(t = 8, r = 6, b = 6, l = 6)
  )

print(p_C)
save_svg(p_C, out_C, width = 4.2, height = 4.2)

# ==============================================================================
# PANEL D — GO enrichment (wide barplot; GO terms on x, angled; colored by module)
# ==============================================================================
enr <- read_tsv(enr_path, show_col_types = FALSE) %>%
  mutate(
    ModuleLabel = factor(as.character(ModuleLabel), levels = MOD_LEVELS),
    GO_Term     = as.character(GO_Term),
    GO_FDR      = as.numeric(GO_FDR),
    GO_neglog10 = as.numeric(GO_neglog10),
    GO_Count    = as.integer(GO_Count)
  ) %>%
  filter(!is.na(GO_Term), GO_Term != "")

TOP_GO_PER_MODULE <- 5

# Select top GO per module (rank by FDR, then by count as tie-breaker)
enr_top <- enr %>%
  group_by(ModuleLabel) %>%
  arrange(GO_FDR, desc(GO_Count)) %>%
  slice_head(n = TOP_GO_PER_MODULE) %>%
  ungroup()

# IMPORTANT: GO terms can repeat across modules; create a unique plotting key
# while still *displaying* GO_Term on the axis.
enr_top <- enr_top %>%
  mutate(
    GO_Key  = paste0(as.character(ModuleLabel), "||", GO_Term),
    GO_Key  = factor(GO_Key, levels = GO_Key)  # temporary; will reset after ordering
  )

# Order x by module then by significance (grouped blocks)
enr_top <- enr_top %>%
  group_by(ModuleLabel) %>%
  arrange(desc(GO_neglog10), .by_group = TRUE) %>%
  mutate(GO_Key = factor(GO_Key, levels = GO_Key)) %>%
  ungroup()

# Vertical separators between module blocks on x (based on GO_Key positions)
bounds <- enr_top %>%
  mutate(xi = as.numeric(GO_Key)) %>%
  group_by(ModuleLabel) %>%
  summarise(xmax = max(xi), .groups = "drop") %>%
  arrange(xmax)

vlines <- bounds$xmax + 0.5
vlines <- vlines[-length(vlines)]

# Short module summaries for legend
mod_summ <- tibble::tibble(
  ModuleLabel = factor(MOD_LEVELS, levels = MOD_LEVELS),
  Summary = c(
    "Protein turnover",
    "Endocytosis & immune interaction",
    "Adaptive immunity",
    "Translation & humoral immunity",
    "ECM remodeling & tissue repair"
  )
)

legend_labs <- setNames(
  paste0(levels(mod_summ$ModuleLabel), ": ", mod_summ$Summary),
  levels(mod_summ$ModuleLabel)
)

p_D <- ggplot(enr_top, aes(x = GO_Key, y = GO_neglog10, fill = ModuleLabel)) +
  geom_col(width = 0.82, color = "grey25", linewidth = 0.25) +
  geom_vline(xintercept = vlines, color = "grey82", linewidth = 0.4) +
  scale_fill_manual(
    values = MOD_COLS,
    breaks = MOD_LEVELS,
    labels = legend_labs,
    name = NULL
  ) +
  # show GO_Term labels on x even though the mapped x is GO_Key
  scale_x_discrete(labels = function(x) sub("^.*\\|\\|", "", x)) +
  labs(
    x = NULL,
    y = expression(paste("-log"[10], "(FDR)"))
  ) +
  theme_qomics() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.35),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 9),
    legend.position = "top",
    plot.margin = margin(t = 14, r = 18, b = 22, l = 18)
  )

print(p_D)
save_svg(p_D, out_D, width = 12.6, height = 3.5)

# ------------------------------------------------------------------------------
# Done
# ------------------------------------------------------------------------------
cat(
  "\n[OK] Figure 6 panels saved in:\n",
  " - ", out_B, "\n",
  " - ", out_C, "\n",
  " - ", out_D, "\n",
  sep = ""
)
