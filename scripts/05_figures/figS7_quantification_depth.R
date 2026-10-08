# ==============================================================================
# Project: Preterm EV Multi-Omics
# Figure:  Fig. S7 – Post-filter quantification depth per sample (proteins, lipids)
# Script:  figS7_quantification_depth.R
# Author:  Nicole Grinovero
#
# Inputs:
#   - data/metadata/sample_metadata.txt
#   - data/proteomics/processed/proteins_filtered.txt  (after fill-rate filter, before imputation)
#   - data/lipidomics/processed/lipids_cleaned.txt     (after fill-rate filter, before imputation)
#
# Output:
#   - outputs/figures/figS7/figS7_quantification_depth.svg
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# ---- Config
tp_levels <- c("T0","T1","T2","T3","T4")

COL_PROT  <- "#8ab6c9"  
COL_LIP   <- "#dfb690"  

# Load metadata
meta <- read_delim(here("data", "metadata", "sample_metadata.txt"), delim = "\t", show_col_types = FALSE) %>%
  mutate(timepoint = factor(timepoint, levels = tp_levels, ordered = TRUE))

# Count non-missing features per sample
count_features <- function(path, meta, name) {
  df <- read_delim(path, delim = "\t", show_col_types = FALSE)
  sample_cols <- intersect(meta$sample_id, names(df))
  counts <- colSums(!is.na(df[, sample_cols]))
  
  tibble(sample_id = names(counts), n_quant = counts) %>%
    left_join(meta %>% select(sample_id, timepoint), by = "sample_id") %>%
    mutate(layer = name)
}

# Load post-filter feature tables
prot_df <- count_features(here("data", "proteomics", "processed", "proteins_filtered.txt"), meta, "Proteins")
lip_df  <- count_features(here("data", "lipidomics", "processed", "lipids_cleaned.txt"), meta, "Lipids")

# --- Combine layers and order samples
barplot_data <- bind_rows(prot_df, lip_df) %>%
  filter(!is.na(timepoint)) %>%
  arrange(timepoint, sample_id) %>%
  mutate(Sample = factor(sample_id, levels = unique(sample_id)))

# Mean per layer (for labels)
avg_stats <- barplot_data %>%
  group_by(layer) %>%
  summarise(y_mean = mean(n_quant), label = paste0("Mean: ", round(y_mean, 0)), .groups = "drop")

# --- Timepoint separators
sep_lines <- barplot_data %>%
  group_by(layer, timepoint) %>%
  summarise(last_idx = max(as.integer(Sample)), .groups = "drop") %>%
  group_by(layer) %>%
  arrange(timepoint) %>%
  filter(timepoint != last(timepoint)) %>% 
  mutate(x_sep = last_idx + 0.5)

# Positions of T0-T4 labels
tp_labels <- barplot_data %>%
  group_by(layer, timepoint) %>%
  summarise(mid_idx = mean(as.integer(Sample)), .groups = "drop")

# ---- Plot
p <- ggplot(barplot_data, aes(x = Sample, y = n_quant, fill = layer)) +
  geom_bar(stat = "identity", width = 0.85) +
  # Mean line
  geom_hline(data = avg_stats, aes(yintercept = y_mean), linetype = "dashed", alpha = 0.5) +
  # Mean label
  geom_text(data = avg_stats, aes(x = 1, y = y_mean, label = label), 
            hjust = 0, vjust = -1, size = 3.5, fontface = "italic", inherit.aes = FALSE) +
  # Timepoint separators
  geom_vline(data = sep_lines, aes(xintercept = x_sep), color = "grey50", linewidth = 0.7) +
  # Timepoint labels
  geom_text(data = tp_labels, aes(x = mid_idx, y = 0, label = timepoint),
            vjust = 2.5, size = 3.5, color = "black", inherit.aes = FALSE) +
  scale_fill_manual(values = c("Proteins" = COL_PROT, "Lipids" = COL_LIP)) +
  facet_wrap(~ layer, scales = "free", ncol = 2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "none",
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold", size = 11),
    plot.margin = margin(t = 10, r = 10, b = 30, l = 10)
  ) +
  labs(x = NULL, y = "Number of non-missing features", title = "Post-filter quantified features per sample") +
  coord_cartesian(clip = "off")

p

dir.create(here("outputs", "figures", "figS7"), recursive = TRUE, showWarnings = FALSE)
ggsave(here("outputs", "figures", "figS7", "figS7_quantification_depth.svg"), plot = p, width = 6.5, height = 3.5)
