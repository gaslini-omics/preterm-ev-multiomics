# ==============================================================================
# Project: Preterm EV Multi-Omics – Proteomics
# Figure:  Fig. S8 – Longitudinal behavior of plasma membrane and endosomal proteins
# Script:  figS8_membrane_endosomal_proteins.R
# Author:  Nicole Grinovero
#
# Aim:
#   Boxplots across T0-T4 of representative transmembrane proteins of the plasma
#   membrane / endosomal compartments (TFRC, LAMP2, ITGB1, CD44) and cytosolic
#   membrane-binding proteins involved in endosomal sorting and vesicle
#   formation (SDCBP, RAB7A, EHD1, SNX2), with their LMM statistics.
#   Inclusion check: significant global time effect (FDR_global < 0.05) and
#   <= 35% imputed values at every timepoint.
#
# Inputs:
#   - data/proteomics/processed/proteins_cleaned.txt       (pre-imputation, for the imputation check)
#   - data/proteomics/processed/proteins_norm_imputed.txt
#   - data/metadata/sample_metadata.txt
#   - outputs/tables/preterm/proteomics/proteins_LMM_time_results.txt
#
# Outputs:
#   - outputs/figures/figS8/figS8_membrane_endosomal_proteins.{svg,png}
#   - outputs/tables/preterm/proteomics/figS8_membrane_endosomal_stats.tsv
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(svglite)
})

# ------------------------------------------------------------------------------
# 0) Paths and settings
# ------------------------------------------------------------------------------
clean_path <- here("data", "proteomics", "processed", "proteins_cleaned.txt")
expr_path  <- here("data", "proteomics", "processed", "proteins_norm_imputed.txt")
meta_path  <- here("data", "metadata", "sample_metadata.txt")
lmm_path   <- here("outputs", "tables", "preterm", "proteomics", "proteins_LMM_time_results.txt")

fig_dir <- here("outputs", "figures", "figS8")
tab_dir <- here("outputs", "tables", "preterm", "proteomics")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

tp_levels   <- c("T0", "T1", "T2", "T3", "T4")
MAX_IMPUTED <- 0.35

panel <- tribble(
  ~Gene,   ~Category,                                    ~Note,
  "TFRC",  "Plasma membrane / endosomal, transmembrane", "transferrin receptor (CD71)",
  "LAMP2", "Plasma membrane / endosomal, transmembrane", "lysosomal/late endosomal membrane",
  "ITGB1", "Plasma membrane / endosomal, transmembrane", "integrin beta-1 (CD29)",
  "CD44",  "Plasma membrane / endosomal, transmembrane", "hyaluronan receptor",
  "SDCBP", "Cytosolic, membrane-binding",                "syntenin-1, ESCRT-accessory",
  "RAB7A", "Cytosolic, membrane-binding",                "late endosomal GTPase",
  "EHD1",  "Cytosolic, membrane-binding",                "endosomal recycling, membrane tubulation",
  "SNX2",  "Cytosolic, membrane-binding",                "retromer, endosomal sorting"
)
cat_levels <- unique(panel$Category)
cat_cols   <- setNames(c("#2F7F7B", "#B5895A"), cat_levels)

# ------------------------------------------------------------------------------
# 1) Load data
# ------------------------------------------------------------------------------
meta <- read_delim(meta_path, delim = "\t", show_col_types = FALSE) %>%
  mutate(timepoint = factor(timepoint, levels = tp_levels))
expr  <- read_delim(expr_path,  delim = "\t", show_col_types = FALSE)
clean <- read_delim(clean_path, delim = "\t", show_col_types = FALSE)
lmm   <- read_delim(lmm_path,   delim = "\t", show_col_types = FALSE)

samples <- intersect(meta$sample_id, names(expr))

to_long <- function(df) {
  df %>%
    filter(GeneName %in% panel$Gene) %>%
    select(GeneName, all_of(samples)) %>%
    pivot_longer(-GeneName, names_to = "sample_id", values_to = "value") %>%
    left_join(meta %>% select(sample_id, timepoint), by = "sample_id")
}

# ------------------------------------------------------------------------------
# 2) Inclusion checks: imputed fraction per timepoint and LMM statistics
# ------------------------------------------------------------------------------
imputed <- to_long(clean) %>%
  group_by(GeneName, timepoint) %>%
  summarise(frac_imputed = mean(is.na(value)), .groups = "drop")

imp_max <- imputed %>%
  group_by(GeneName) %>%
  summarise(max_frac_imputed = max(frac_imputed), .groups = "drop")

stats <- panel %>%
  left_join(lmm %>% select(GeneName, beta_T4_vs_T0, FDR_T4_T0, FDR_global),
            by = c("Gene" = "GeneName")) %>%
  left_join(imp_max, by = c("Gene" = "GeneName")) %>%
  mutate(quantified = !is.na(beta_T4_vs_T0))

print(as.data.frame(stats))
stopifnot(all(stats$FDR_global < 0.05), all(stats$max_frac_imputed <= MAX_IMPUTED))

write_tsv(stats, file.path(tab_dir, "figS8_membrane_endosomal_stats.tsv"))

# ------------------------------------------------------------------------------
# 3) Plot
# ------------------------------------------------------------------------------
plot_df <- to_long(expr) %>%
  left_join(panel, by = c("GeneName" = "Gene")) %>%
  mutate(GeneName = factor(GeneName, levels = panel$Gene),
         Category = factor(Category, levels = cat_levels))

p <- ggplot(plot_df, aes(timepoint, value)) +
  geom_boxplot(aes(fill = Category), colour = "#595959", alpha = 0.6,
               width = 0.6, outlier.shape = NA, linewidth = 0.6) +
  geom_point(colour = "#B3B3B3", alpha = 0.6, size = 1.2,
             position = position_jitter(width = 0.12, height = 0, seed = 1)) +
  stat_summary(fun = median, geom = "point", shape = 21, fill = "white",
               colour = "#333333", size = 1.6, stroke = 0.4) +
  facet_wrap(~ GeneName, nrow = 2, scales = "free_y") +
  scale_fill_manual(values = cat_cols, name = NULL) +
  labs(x = NULL, y = expression(log[2] * "(abundance)")) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid       = element_blank(),
    strip.background = element_rect(fill = "#D9D9D9"),
    strip.text       = element_text(size = 11),
    legend.position  = "bottom",
    text             = element_text(family = "Arial")
  )

ggsave(file.path(fig_dir, "figS8_membrane_endosomal_proteins.svg"), p,
       width = 9, height = 5, device = svglite::svglite)
ggsave(file.path(fig_dir, "figS8_membrane_endosomal_proteins.png"), p,
       width = 9, height = 5, dpi = 600, bg = "white")

message("Fig. S8 written to: ", fig_dir)
