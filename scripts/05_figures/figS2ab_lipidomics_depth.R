# ==============================================================================
# Fig. S2a,b — Lipidomics method validation (cutoff rationale + depth without cutoff)
#   (Fig. S2c: fig1_method_validation.R; Fig. S2d: fig1d_figS2d_ev_signature.R)
# Project:   Mag-Net EV Multi-Omics — Technical validation (tech reps)
# Author:    Nicole Grinovero
# Date:      Feb 2026
#
# Purpose:
#   a) Lipid intensity distributions (POS/NEG) with cutoff line at 1e5 (log10 scale)
#   b) Lipidomics depth WITHOUT intensity cutoff (Intensity > 0), POS/NEG
#
# Inputs:
#   - data/techrep/lipidomics/processed/techrep_lipids_cleaned_C18_pos.txt
#   - data/techrep/lipidomics/processed/techrep_lipids_cleaned_C18_neg.txt
#   - data/techrep/lipidomics/sample_meta_C18_pos.txt
#   - data/techrep/lipidomics/sample_meta_C18_neg.txt
#
# Outputs:
#   - outputs/figures/figS2/
#     - figS2a_lipid_intensity_distributions_cutoff.(svg/png)
#     - figS2b_lipidomics_depth_no_cutoff.(svg/png)
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(stringr)
  library(patchwork)
  library(svglite)
  library(scales)
})

# ------------------------------------------------------------------------------
# 0) PATHS
# ------------------------------------------------------------------------------
in_pos   <- here("data","techrep","lipidomics","processed","techrep_lipids_cleaned_C18_pos.txt")
in_neg   <- here("data","techrep","lipidomics","processed","techrep_lipids_cleaned_C18_neg.txt")
meta_pos <- here("data","techrep","lipidomics","sample_meta_C18_pos.txt")
meta_neg <- here("data","techrep","lipidomics","sample_meta_C18_neg.txt")

out_dir  <- here("outputs","figures","figS2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1) SETTINGS
# ------------------------------------------------------------------------------
cutoff_intensity <- 1e5
lip_presence_threshold_no_cutoff <- 0  # Intensity > 0 for "no cutoff" depth

# Fixed logical condition order (used as a fallback + for stable labeling)
COND_ORDER <- c(
  "QC",
  "Mag-Net BAMM 10µL","Mag-Net BAMM 5µL",
  "Mag-Net IPA 10µL","Mag-Net IPA 5µL",
  "Mag-Net MeOH 10µL","Mag-Net MeOH 5µL",
  "Mag-Net standard 10µL","Mag-Net standard 5µL",
  "Neat"
)

# ------------------------------------------------------------------------------
# 2) PALETTE + THEME (same as fig1_method_validation.R)
# ------------------------------------------------------------------------------
cond_colors <- c(
  "Neat"                  = "whitesmoke",
  "Mag-Net standard 5µL"  = "#B8B8B8",
  "Mag-Net standard 10µL" = "#9E9E9E",
  "Mag-Net MeOH 5µL"      = "#F3E8C9",
  "Mag-Net MeOH 10µL"     = "#E6D3A3",
  "Mag-Net BAMM 5µL"      = "#C7DDE8",
  "Mag-Net BAMM 10µL"     = "#9BBFD4",
  "Mag-Net IPA 5µL"       = "#E1D7E6",
  "Mag-Net IPA 10µL"      = "#C3B2D1",
  "QC"                    = "#EFEFEF",
  "BLK"                   = "#FFFFFF"
)

BASE_SIZE <- 12

theme_qomics <- function() {
  theme_bw(base_size = BASE_SIZE) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color = "grey10"),
      axis.title = element_text(color = "grey10"),
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey92", color = NA),
      plot.title = element_text(size = 14),
      plot.margin = margin(6, 26, 6, 6),  # extra right margin (numbers inside bars)
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 9)
    )
}

save_svg_png <- function(p, prefix, w, h, dpi = 320) {
  ggsave(paste0(prefix, ".svg"), p, width = w, height = h, device = svglite::svglite)
  ggsave(paste0(prefix, ".png"), p, width = w, height = h, dpi = dpi, bg = "white")
}

# ------------------------------------------------------------------------------
# 3) HELPERS
# ------------------------------------------------------------------------------
norm_ul <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("μL", "µL") %>%
    str_replace_all("(?i)uL", "µL") %>%
    str_replace_all("\\s+", " ") %>%
    str_squish()
}

is_detected_no_cutoff <- function(x, thr = 0) {
  x <- tidyr::replace_na(x, 0)
  is.finite(x) & x > thr
}

parse_condition_tbl <- function(df, condition_col = "Condition") {
  cond <- df[[condition_col]] %>% as.character() %>% norm_ul()
  
  vol <- str_extract(cond, "(?i)(5|10)\\s*µL") %>%
    str_replace_all("[^0-9]", "") %>%
    suppressWarnings() %>%
    as.integer()
  
  base <- cond %>%
    str_replace_all("(?i)\\s*(5|10)\\s*µL\\s*$", "") %>%
    str_squish()
  
  df %>%
    mutate(
      Condition = cond,
      BaseCond  = base,
      Vol_uL    = vol
    )
}

# Facet-specific ordering (stable):
# - rank BaseCond within each Polarity by BaseMean = mean(depth across 5+10)
# - within BaseCond: 10µL before 5µL
# - QC always first if present
make_facet_levels <- function(depth_sum) {
  base_rank <- depth_sum %>%
    group_by(Polarity, BaseCond) %>%
    summarise(BaseMean = mean(Mean, na.rm = TRUE), .groups = "drop") %>%
    group_by(Polarity) %>%
    arrange(desc(BaseMean), .by_group = TRUE) %>%
    summarise(base_ord = list(BaseCond), .groups = "drop")
  
  cond_levels_by_pol <- depth_sum %>%
    left_join(base_rank, by = "Polarity") %>%
    group_by(Polarity) %>%
    summarise(
      cond_levels = list({
        df <- cur_data() %>%
          mutate(
            BaseCond = factor(BaseCond, levels = base_ord[[1]]),
            VolSort  = if_else(is.na(Vol_uL), 999L, -Vol_uL) # 10 before 5
          ) %>%
          arrange(BaseCond, VolSort)
        
        lev <- df$Condition
        lev <- c("QC", setdiff(lev, "QC"))
        unique(lev)
      }),
      .groups = "drop"
    )
  
  cond_levels_by_pol
}

apply_panel_levels <- function(df, cond_levels_by_pol) {
  get_levels <- function(pol) cond_levels_by_pol$cond_levels[[match(pol, cond_levels_by_pol$Polarity)]]
  
  df %>%
    mutate(
      Condition = norm_ul(Condition),
      CondPanel = paste0(as.character(Polarity), "__", Condition)
    ) %>%
    group_by(Polarity) %>%
    group_modify(~{
      pol <- as.character(.y$Polarity)
      lev <- get_levels(pol)
      lev_panel <- paste0(pol, "__", lev)
      .x %>% mutate(CondPanel = factor(CondPanel, levels = unique(lev_panel)))
    }) %>%
    ungroup()
}

# ------------------------------------------------------------------------------
# 4) LOAD META
# ------------------------------------------------------------------------------
metaP <- readr::read_tsv(meta_pos, show_col_types = FALSE) %>%
  mutate(
    Condition = norm_ul(Condition),
    Sample_ID = as.character(Sample_ID)
  )
metaN <- readr::read_tsv(meta_neg, show_col_types = FALSE) %>%
  mutate(
    Condition = norm_ul(Condition),
    Sample_ID = as.character(Sample_ID)
  )

# ------------------------------------------------------------------------------
# 5) LOAD MATRICES + LONG FORMAT
# ------------------------------------------------------------------------------
matP <- readr::read_tsv(in_pos, show_col_types = FALSE) %>% mutate(`Alignment ID` = as.character(`Alignment ID`))
matN <- readr::read_tsv(in_neg, show_col_types = FALSE) %>% mutate(`Alignment ID` = as.character(`Alignment ID`))

sampP <- intersect(colnames(matP), metaP$Sample_ID)
sampN <- intersect(colnames(matN), metaN$Sample_ID)

longP <- matP %>%
  mutate(.feature_id = row_number()) %>%
  select(all_of(sampP), .feature_id) %>%
  pivot_longer(cols = all_of(sampP), names_to = "Sample_ID", values_to = "Intensity") %>%
  left_join(metaP, by = "Sample_ID") %>%
  mutate(Polarity = "Positive polarity")

longN <- matN %>%
  mutate(.feature_id = row_number()) %>%
  select(all_of(sampN), .feature_id) %>%
  pivot_longer(cols = all_of(sampN), names_to = "Sample_ID", values_to = "Intensity") %>%
  left_join(metaN, by = "Sample_ID") %>%
  mutate(Polarity = "Negative polarity")

lip_all <- bind_rows(longP, longN) %>%
  mutate(
    Condition = norm_ul(Condition),
    Polarity  = factor(Polarity, levels = c("Negative polarity", "Positive polarity"))
  )

# ------------------------------------------------------------------------------
# 6) PANEL a — Intensity distributions + cutoff line
# ------------------------------------------------------------------------------
pS2a <- lip_all %>%
  filter(!is.na(Intensity), is.finite(Intensity), Intensity > 0) %>%
  ggplot(aes(x = log10(Intensity))) +
  geom_histogram(bins = 90, fill = "grey70", color = "grey70") +
  geom_vline(xintercept = log10(cutoff_intensity), linetype = 2, linewidth = 0.9, color = "grey15") +
  facet_wrap(~Polarity, nrow = 1, scales = "free_y") +
  labs(
    title = "Lipid intensity distributions",
    x = expression(log[10]~"(Intensity)"),
    y = "Count"
  ) +
  theme_qomics() +
  theme(
    axis.title.y = element_text(size = 12),
    plot.margin = margin(10, 20, 10, 12)
  )

pS2a
save_svg_png(pS2a, file.path(out_dir, "figS2a_lipid_intensity_distributions_cutoff"), w = 5, h = 3)

# ------------------------------------------------------------------------------
# 7) PANEL b — Lipidomics depth WITHOUT cutoff (Intensity > 0)
# ------------------------------------------------------------------------------
depth_by_sample_no <- lip_all %>%
  filter(!is.na(Condition), Condition != "BLK") %>%
  group_by(Polarity, Sample_ID, Condition) %>%
  summarise(
    n_detected = sum(is_detected_no_cutoff(Intensity, thr = lip_presence_threshold_no_cutoff)),
    .groups = "drop"
  ) %>%
  parse_condition_tbl()

depth_sum_no <- depth_by_sample_no %>%
  group_by(Polarity, Condition, BaseCond, Vol_uL) %>%
  summarise(
    Mean = mean(n_detected, na.rm = TRUE),
    SD   = sd(n_detected, na.rm = TRUE),
    .groups = "drop"
  )

# facet-specific x ordering
cond_levels_by_pol <- make_facet_levels(depth_sum_no)
depth_sum_no       <- apply_panel_levels(depth_sum_no, cond_levels_by_pol)
depth_by_sample_no <- apply_panel_levels(depth_by_sample_no, cond_levels_by_pol)

pS2b <- ggplot(depth_sum_no, aes(x = CondPanel, y = Mean, fill = Condition)) +
  geom_col(width = 0.78, alpha = 0.95, color = "grey25", linewidth = 0.25) +
  geom_errorbar(
    aes(ymin = pmax(0, Mean - SD), ymax = Mean + SD),
    width = 0.22, linewidth = 0.30, color = "grey45"
  ) +
  geom_point(
    data = depth_by_sample_no,
    aes(x = CondPanel, y = n_detected),
    position = position_jitter(width = 0.12, height = 0),
    size = 2.0, alpha = 0.30, color = "grey10",
    inherit.aes = FALSE
  ) +
  geom_text(
    aes(label = as.integer(round(Mean))),
    vjust = -0.55, size = 3.4, color = "grey20"
  ) +
  facet_wrap(~Polarity, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = cond_colors, drop = FALSE) +
  scale_x_discrete(labels = function(x) sub("^.*__", "", x)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "Lipidomics depth without intensity cutoff",
    x = NULL,
    y = "Detected lipid features"
  ) +
  theme_qomics() +
  theme(plot.margin = margin(t = 15, r = 20, b = 20, l = 20),
        legend.position = "none",
        axis.text.x = element_text(angle = 50, hjust = 1)) +
  coord_cartesian(clip = "off")

pS2b
save_svg_png(pS2b, file.path(out_dir, "figS2b_lipidomics_depth_no_cutoff"), w = 5, h = 4)

# ------------------------------------------------------------------------------
# 8) COMBINE (S2a over S2b)
# ------------------------------------------------------------------------------
FigS2 <- pS2a / pS2b + plot_layout(heights = c(1, 1.05))
FigS2
save_svg_png(FigS2, file.path(out_dir, "figS2ab_combined"), w = 10.2, h = 8.4)

# ------------------------------------------------------------------------------
# 9) SOURCE DATA
# ------------------------------------------------------------------------------
write_tsv(
  depth_by_sample_no %>%
    mutate(Condition = norm_ul(Condition)) %>%
    select(Polarity, Sample_ID, Condition, BaseCond, Vol_uL, n_detected),
  file.path(out_dir, "figS2b_source_data.tsv")
)

message("Done. Outputs in: ", out_dir)
