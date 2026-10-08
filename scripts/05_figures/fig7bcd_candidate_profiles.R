# ==============================================================================
# Project: Preterm EV Multi-Omics – Figure 7b–d (candidate profiles)
# Script:  fig7bcd_candidate_profiles.R
# Author:  Nicole Grinovero
# Aim:     Candidate-profile panels (leading-edge features) for Fig. 7b–d.
#          Mean trajectories ± SE by clinical group.
#
# Inputs:
#   - data/metadata/sample_metadata.txt
#   - data/proteomics/processed/proteins_norm_imputed.txt
#   - data/lipidomics/processed/lipids_norm_imputed.txt
#   - outputs/tables/preterm/integration/integration_modules_membership_annotated.tsv
#
# Output:
#   - outputs/figures/fig7_clinical_associations/fig7bcd_candidates.svg
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# 1) Setup
# ------------------------------------------------------------------------------
fig_dir <- here("outputs", "figures", "fig7_clinical_associations")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Group colours
COLOR_M1_HEM    <- c("Control" = "grey70", "Brain lesion (hemorrhagic)" = "#1F4E79")
COLOR_M3_BRAIN  <- c("Control" = "grey70", "Brain lesion (any)"         = "#6B8EC6")
COLOR_M3_SEPSIS <- c("Control" = "grey70", "Sepsis"                    = "#C26A3D" )

STRIP_M1 <- "#C8D7C1"
STRIP_M3 <- "#C9BEDD"

# Paths
meta_path      <- here("data", "metadata", "sample_metadata.txt")
prot_expr_path <- here("data", "proteomics", "processed", "proteins_norm_imputed.txt")
lip_expr_path  <- here("data", "lipidomics",  "processed", "lipids_norm_imputed.txt")
memb_path      <- here("outputs", "tables", "preterm", "integration",
                       "integration_modules_membership_annotated.tsv")

# ------------------------------------------------------------------------------
# 2) Load data
# ------------------------------------------------------------------------------
meta <- read_delim(meta_path, delim = "\t", show_col_types = FALSE) %>%
  mutate(timepoint = factor(timepoint, levels = c("T0","T1","T2","T3","T4"), ordered = TRUE))

annot <- read_delim(memb_path, delim = "\t", show_col_types = FALSE)
if (!"Feature" %in% names(annot)) colnames(annot)[1] <- "Feature"

annot <- annot %>%
  select(Feature, GeneName, Type) %>%
  distinct(Feature, .keep_all = TRUE)

load_mat <- function(path, id_col) {
  read_delim(path, delim = "\t", show_col_types = FALSE) %>%
    select(all_of(id_col), where(is.numeric)) %>%
    distinct(.data[[id_col]], .keep_all = TRUE) %>%
    column_to_rownames(id_col) %>%
    as.matrix()
}

mat_prot <- load_mat(prot_expr_path, "UniProtID")
mat_lip  <- load_mat(lip_expr_path,  "LipidName")

find_id <- function(gene_symbol) {
  res <- annot %>% filter(GeneName == gene_symbol) %>% pull(Feature)
  if (length(res) == 0) NA_character_ else res[1]
}

# ------------------------------------------------------------------------------
# 3) Features shown (leading-edge features of the significant module–outcome
#    associations; see drivers_top2_<layer>_<outcome>.tsv from 04_integration/03)
# ------------------------------------------------------------------------------
# Panel B — M1 proteins (hemorrhagic lesion)
ids_m1 <- c(find_id("C4BPA"), find_id("C4BPB"), find_id("HBM"))
ids_m1 <- ids_m1[!is.na(ids_m1)]

# Panel C — M3 proteins (any brain lesion)
ids_m3_prot <- c(find_id("GCA"), find_id("NIBAN1"), find_id("PSMB8"))
ids_m3_prot <- ids_m3_prot[!is.na(ids_m3_prot)]

# Panel D — M3 lipids (sepsis)
# Must match rownames(mat_lip)
ids_m3_lip <- c("PC 33:1", "PC 32:1")

# ------------------------------------------------------------------------------
# 4) Plot function
# ------------------------------------------------------------------------------
plot_clean <- function(target_ids, mat, meta_df,
                       outcome_col, outcome_label,
                       color_map,
                       title, strip_col,
                       ncol_facets = 3) {
  
  hits <- intersect(target_ids, rownames(mat))
  if (length(hits) == 0) return(NULL)
  
  d <- t(mat[hits, , drop = FALSE]) %>%
    as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    pivot_longer(-sample_id, names_to = "Feature", values_to = "Expression") %>%
    inner_join(meta_df, by = "sample_id") %>%
    left_join(annot, by = "Feature") %>%
    mutate(
      DisplayName = ifelse(!is.na(GeneName) & nzchar(GeneName), GeneName, Feature),
      Group = ifelse(.data[[outcome_col]] == 1, outcome_label, "Control")
    ) %>%
    filter(!is.na(timepoint), !is.na(Expression), !is.na(Group))
  
  # enforce facet order based on target_ids
  display_levels <- tibble(Feature = target_ids) %>%
    left_join(annot, by = "Feature") %>%
    mutate(DisplayName = ifelse(!is.na(GeneName) & nzchar(GeneName), GeneName, Feature)) %>%
    pull(DisplayName) %>%
    unique()
  
  # for lipids: if no annot/gene, keep Feature order
  if (length(display_levels) == 0 || all(is.na(display_levels))) {
    display_levels <- unique(target_ids)
  }
  
  d <- d %>%
    mutate(DisplayName = factor(DisplayName, levels = display_levels))
  
  # counts (for subtitle)
  n_case <- sum(d$Group == outcome_label, na.rm = TRUE)
  n_ctrl <- sum(d$Group == "Control", na.rm = TRUE)
  
  ggplot(d, aes(x = timepoint, y = Expression, color = Group, group = Group)) +
    
    # mean line
    stat_summary(fun = mean, geom = "line", linewidth = 1.8) +
    stat_summary(fun = mean, geom = "point", size = 2.4) +
    
    # mean ± SE
    stat_summary(
      fun.data = mean_se,
      geom = "errorbar",
      width = 0.18,
      linewidth = 0.7
    ) +
    
    facet_wrap(~DisplayName, scales = "free_y", ncol = ncol_facets) +
    
    scale_color_manual(values = color_map) +
    
    labs(
      title = title,
      subtitle = paste0(outcome_label, " n=", n_case, " | Control n=", n_ctrl),
      y = expression("log"[2]*"(intensity)"),
      x = NULL,
      color = NULL
    ) +
    
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "grey30"),
      strip.background = element_rect(fill = strip_col, color = NA),
      strip.text = element_text(face = "bold", size = 10),
      legend.position = "bottom",
      panel.spacing = unit(1.2, "lines")
    )
}

# ------------------------------------------------------------------------------
# 5) Generate panels
# ------------------------------------------------------------------------------
p_m1 <- plot_clean(
  target_ids = ids_m1,
  mat = mat_prot,
  meta_df = meta,
  outcome_col = "hemorrhagic_lesion",
  outcome_label = "Brain lesion (hemorrhagic)",
  color_map = COLOR_M1_HEM,
  title = "M1 proteins — leading candidates (hemorrhagic lesion)",
  strip_col = STRIP_M1,
  ncol_facets = 3
)

p_m3 <- plot_clean(
  target_ids = ids_m3_prot,
  mat = mat_prot,
  meta_df = meta,
  outcome_col = "any_brain_lesion",
  outcome_label = "Brain lesion (any)",
  color_map = COLOR_M3_BRAIN,
  title = "M3 proteins — leading candidates (any brain lesion)",
  strip_col = STRIP_M3,
  ncol_facets = 3
)

p_lip <- plot_clean(
  target_ids = ids_m3_lip,
  mat = mat_lip,
  meta_df = meta,
  outcome_col = "sepsis",
  outcome_label = "Sepsis",
  color_map = COLOR_M3_SEPSIS,
  title = "M3 lipids — leading candidates (sepsis)",
  strip_col = STRIP_M3,
  ncol_facets = 2
)

# Row C: 2 lipid panels + spacer
lip_row <- (p_lip + plot_spacer()) +
  plot_layout(widths = c(1, 0.5))

final_fig7 <- (p_m1 / p_m3 / lip_row) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

final_fig7

# ------------------------------------------------------------------------------
# Export SVG (vector)
# ------------------------------------------------------------------------------
out_svg <- file.path(fig_dir, "fig7bcd_candidates.svg")

ggsave(
  filename = out_svg,
  plot = final_fig7,
  device = svglite::svglite,
  width = 8.5,
  height = 8.5,
  units = "in"
)
