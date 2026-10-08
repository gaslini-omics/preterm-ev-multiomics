# ==============================================================================
# Fig. 1 — Technical validation of the EV multi-omics workflow
# Script panel -> manuscript panel
#   A) Proteomics depth (Peptides / Protein groups)            -> Fig. 1a
#   B) Proteomics reproducibility (CV bins, 100% stacked)      -> Fig. 1b
#   C) Proteomics comparability vs Mag-Net standard (scatter)  -> Fig. 1c
#   D) Lipidomics depth (Intensity > 1e5), POS/NEG             -> Fig. 1e
#   E) Lipid class composition, POS/NEG                        -> Fig. S2c
#   F) Lipidomics reproducibility (CV ridge), POS/NEG          -> Fig. 1f
#   (Fig. 1d is produced by fig1d_figS2d_ev_signature.R)
#
# Inputs: data/techrep/proteomics/ (Spectronaut), data/techrep/lipidomics/
#         (output of 01_techrep_lip_cleaning.R + sample metadata)
# Output: outputs/figures/fig1_technical_validation/
#
# Author: Nicole Grinovero
# Date: Feb 2026
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
# Proteomics
idov_path   <- here("data","techrep","proteomics","plasma_EV_multiomics_techrep_IdentificationsOverview.tsv")
report_path <- here("data","techrep","proteomics","20260206_182115_plasma_EV_multiomics_techrep_Report.tsv")

# Lipidomics
in_pos   <- here("data","techrep","lipidomics","processed","techrep_lipids_cleaned_C18_pos.txt")
in_neg   <- here("data","techrep","lipidomics","processed","techrep_lipids_cleaned_C18_neg.txt")
meta_pos <- here("data","techrep","lipidomics","sample_meta_C18_pos.txt")
meta_neg <- here("data","techrep","lipidomics","sample_meta_C18_neg.txt")

# Output
out_main <- here("outputs","figures","fig1_technical_validation")
dir.create(out_main, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1) PALETTES
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

cv_group_colors <- c(
  "0–10%"  = "#9FC5A8",
  "10–20%" = "#D8E2A8",
  "20–30%" = "#F0D7A8",
  ">30%"   = "#D9A3A3"
)

solv_cols <- c(
  "BAMM" = "#9BBFD4",
  "IPA"  = "#C3B2D1",
  "MeOH" = "#E6D3A3"
)

# ------------------------------------------------------------------------------
# 2) STYLE
# ------------------------------------------------------------------------------
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
    str_replace_all("μL", "µL") %>%          # greek mu -> micro sign
    str_replace_all("(?i)uL", "µL") %>%      # uL -> µL
    str_replace_all("\\s+", " ") %>%         # collapse spaces
    str_squish()
}

cutoff_intensity <- 1e5
is_detected <- function(x) is.finite(x) & !is.na(x) & x > cutoff_intensity

recode_condition <- function(x) {
  case_when(
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
}

# fixed order of solvents (project-driven)
SOLVENT_ORDER <- c("BAMM", "IPA", "MeOH", "standard", "Other")

parse_solvent <- function(cond) {
  cond <- as.character(cond)
  case_when(
    str_detect(cond, "BAMM")     ~ "BAMM",
    str_detect(cond, "IPA")      ~ "IPA",
    str_detect(cond, "MeOH")     ~ "MeOH",
    str_detect(cond, "standard") ~ "standard",
    TRUE                         ~ "Other"
  )
}

parse_volume <- function(cond) {
  cond <- as.character(cond)
  case_when(
    str_detect(cond, "10\\s*µL") ~ "10",
    str_detect(cond, "5\\s*µL")  ~ "5",
    TRUE                         ~ NA_character_
  )
}

# --- PROTEOMICS ORDER: BAMM -> IPA -> MeOH -> standard, 10 -> 5, Neat last ---
build_order_conditions <- function(conds) {
  tb <- tibble(Condition = unique(as.character(conds))) %>%
    mutate(
      Condition = norm_ul(Condition),
      solvent   = parse_solvent(Condition),
      volume    = parse_volume(Condition)
    )
  
  # keep Neat always last (but present if exists)
  neat_present <- any(tb$Condition == "Neat", na.rm = TRUE)
  
  ord <- tb %>%
    filter(Condition != "Neat") %>%
    mutate(
      solvent = factor(solvent, levels = SOLVENT_ORDER),
      volume  = factor(volume,  levels = c("10", "5"))
    ) %>%
    arrange(solvent, volume) %>%
    pull(Condition) %>%
    unique()
  
  if (neat_present) c(ord, "Neat") else ord
}

# ------------------------------------------------------------------------------
# PANEL A — Proteomics depth (IdentificationsOverview.tsv)
#   vertical bars (like lipidomics), fixed solvent order, numbers inside bars
# ------------------------------------------------------------------------------
id_data <- read_tsv(idov_path, show_col_types = FALSE)

id_long <- id_data %>%
  transmute(
    Condition = norm_ul(recode_condition(Condition)),
    Replicate,
    Peptides,
    `Protein groups` = ProteinGroups
  ) %>%
  pivot_longer(cols = c("Peptides", "Protein groups"),
               names_to = "Metric", values_to = "IDs")

ord_levels <- build_order_conditions(id_long$Condition)

id_long <- id_long %>%
  mutate(CondPanel = factor(Condition, levels = ord_levels))

id_sum <- id_long %>%
  group_by(Metric, CondPanel) %>%
  summarise(
    Mean = mean(IDs, na.rm = TRUE),
    SD   = sd(IDs, na.rm = TRUE),
    .groups = "drop"
  )

pA <- ggplot(id_sum, aes(x = CondPanel, y = Mean, fill = CondPanel)) +
  geom_col(width = 0.78, alpha = 0.95, color = "grey25", linewidth = 0.25) +
  geom_errorbar(aes(ymin = pmax(0, Mean - SD), ymax = Mean + SD),
                width = 0.22, linewidth = 0.30, color = "grey45") +
  geom_point(
    data = id_long,
    aes(x = CondPanel, y = IDs),
    position = position_jitter(width = 0.12, height = 0),
    size = 2.0, alpha = 0.30, color = "grey10",
    inherit.aes = FALSE
  ) +
  geom_text(
    aes(label = as.integer(round(Mean, 0))),
    angle = 90, vjust = 0.5, hjust = 1.3,
    size = 3.2, color = "grey10"
  ) +
  facet_wrap(~Metric, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = cond_colors, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Proteomics depth", x = NULL, y = "Mean IDs") +
  theme_qomics() +
  theme(plot.margin = margin(t = 15, r = 20, b = 20, l = 20),
        axis.text.x = element_text(angle = 50, hjust = 1),
        legend.position = "none") +
  coord_cartesian(clip = "off")

pA

save_svg_png(pA, file.path(out_main, "fig1a_proteomics_depth"), w = 6.3, h = 4)

# ------------------------------------------------------------------------------
# PANEL B — Proteomics CV distribution (Report.tsv) — 100% stacked
# ------------------------------------------------------------------------------
rep <- read_delim(report_path, delim = "\t", show_col_types = FALSE)
stopifnot(all(c("PG.ProteinGroups","PG.Genes") %in% colnames(rep)))

sample_cols <- grep("\\.PG\\.Quantity$", colnames(rep), value = TRUE)
stopifnot(length(sample_cols) > 0)

sample_meta <- tibble(sample_col = sample_cols) %>%
  mutate(
    sample_clean  = str_remove(sample_col, "^\\[\\d+\\]\\s+"),
    replicate     = str_extract(sample_clean, "C\\d+"),
    condition_raw = str_remove(sample_clean, "_C\\d+_\\d+\\.raw\\.PG\\.Quantity$"),
    Condition     = norm_ul(recode_condition(condition_raw)),
    Condition     = factor(Condition, levels = ord_levels)
  )

expr_long <- rep %>%
  select(PG.ProteinGroups, PG.Genes, all_of(sample_cols)) %>%
  pivot_longer(cols = all_of(sample_cols), names_to = "sample_col", values_to = "Intensity") %>%
  left_join(sample_meta, by = "sample_col") %>%
  mutate(
    Intensity = as.numeric(Intensity),
    Present   = is.finite(Intensity) & !is.na(Intensity)
  )

cv_tbl <- expr_long %>%
  group_by(Condition, PG.ProteinGroups) %>%
  summarise(
    n_obs = sum(Present),
    meanI = mean(Intensity, na.rm = TRUE),
    sdI   = sd(Intensity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(CV = if_else(n_obs >= 2 & is.finite(meanI) & meanI > 0, 100 * sdI / meanI, NA_real_)) %>%
  filter(is.finite(CV))

cv_binned <- cv_tbl %>%
  mutate(
    CV_group = case_when(
      CV <= 10 ~ "0–10%",
      CV <= 20 ~ "10–20%",
      CV <= 30 ~ "20–30%",
      TRUE     ~ ">30%"
    ),
    CV_group = factor(CV_group, levels = c("0–10%", "10–20%", "20–30%", ">30%")),
    # enforce fixed condition order (same as Panel A)
    Condition = factor(as.character(Condition), levels = ord_levels)
  ) %>%
  count(Condition, CV_group, name = "n") %>%
  group_by(Condition) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

pB <- ggplot(cv_binned, aes(x = Condition, y = prop, fill = CV_group)) +
  geom_col(width = 0.86, color = "white", linewidth = 0.5) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_manual(values = cv_group_colors, name = "CV group") +
  labs(
    title = "Proteomics reproducibility (CV distribution)",
    x = NULL, y = "Proportion of protein groups (%)"
  ) +
  theme_qomics() +
  theme(
    plot.margin = margin(t = 15, r = 20, b = 20, l = 20),
    axis.text.x = element_text(angle = 50, hjust = 1),
    legend.position = "right"
  )

pB
save_svg_png(pB, file.path(out_main, "fig1b_proteomics_cv"), w = 4.8, h = 4)

# ------------------------------------------------------------------------------
# PANEL C — Proteomics comparability vs Mag-Net standard (scatter)
# ------------------------------------------------------------------------------
ref_10   <- "Mag-Net standard 10µL"
ref_5    <- "Mag-Net standard 5µL"
conds_10 <- c("Mag-Net BAMM 10µL", "Mag-Net IPA 10µL", "Mag-Net MeOH 10µL")
conds_5  <- c("Mag-Net BAMM 5µL",  "Mag-Net IPA 5µL",  "Mag-Net MeOH 5µL")

# subset from expr_long already built
expr_C <- expr_long %>%
  mutate(Condition = as.character(Condition)) %>%
  filter(Condition %in% c(ref_10, ref_5, conds_10, conds_5)) %>%
  filter(Condition != "Neat") %>%  # exclude neat in comparability
  mutate(
    Gene = str_squish(as.character(PG.Genes)),
    Log2 = log2(Intensity + 1)
  )

# one protein group per gene (stable)
gene_pick <- expr_C %>%
  group_by(Gene, PG.ProteinGroups) %>%
  summarise(med = median(Intensity, na.rm = TRUE), .groups = "drop") %>%
  group_by(Gene) %>%
  slice_max(med, n = 1, with_ties = FALSE) %>%
  ungroup()

expr_C <- expr_C %>% semi_join(gene_pick, by = c("Gene", "PG.ProteinGroups"))

cond_mean <- expr_C %>%
  group_by(Gene, Condition) %>%
  summarise(Log2 = mean(Log2, na.rm = TRUE), .groups = "drop")

wide <- cond_mean %>% pivot_wider(names_from = Condition, values_from = Log2)

make_scatter <- function(ref_cond, tests, vol_label) {
  wide_small <- wide %>% select(Gene, all_of(c(ref_cond, tests)))
  
  sc <- wide_small %>%
    pivot_longer(cols = all_of(tests), names_to = "TestCondition", values_to = "TestLog2") %>%
    mutate(
      RefLog2 = .data[[ref_cond]],
      Volume  = vol_label
    ) %>%
    filter(is.finite(RefLog2), is.finite(TestLog2)) %>%
    left_join(
      tibble(TestCondition = tests, Solvent = parse_solvent(tests)),
      by = "TestCondition"
    ) %>%
    mutate(Solvent = factor(Solvent, levels = c("BAMM","IPA","MeOH"))) %>%
    select(Gene, Volume, Solvent, RefLog2, TestLog2)
  
  labs <- sc %>%
    group_by(Volume, Solvent) %>%
    summarise(r = cor(RefLog2, TestLog2, use = "pairwise.complete.obs", method = "pearson"),
              .groups = "drop") %>%
    mutate(label = sprintf("r = %.2f", r))
  
  list(sc = sc, labs = labs)
}

d10 <- make_scatter(ref_10, conds_10, "10 µL input")
d5  <- make_scatter(ref_5,  conds_5,  "5 µL input")

plot_df <- bind_rows(d10$sc, d5$sc) %>%
  mutate(Volume = factor(Volume, levels = c("10 µL input", "5 µL input")))

lab_df <- bind_rows(d10$labs, d5$labs) %>%
  mutate(Volume = factor(Volume, levels = c("10 µL input", "5 µL input")))

pC <- ggplot(plot_df, aes(x = RefLog2, y = TestLog2)) +
  geom_point(aes(color = Solvent), alpha = 0.5, size = 1.25) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.7, color = "grey20") +
  facet_grid(Volume ~ Solvent) +
  geom_text(
    data = lab_df,
    aes(x = -Inf, y = Inf, label = label),
    hjust = -0.10, vjust = 1.15, size = 3.4,
    inherit.aes = FALSE
  ) +
  coord_equal() +
  labs(
    title = "Proteomics comparability to standard Mag-Net",
    x = expression("Mag-Net standard (log"[2]*" intensity)"),
    y = expression("Test condition (log"[2]*" intensity)")
  ) +
  theme_qomics() +
  theme(plot.margin = margin(t = 15, r = 20, b = 20, l = 20),
        legend.position = "none") +
  scale_color_manual(values = solv_cols)

pC

save_svg_png(pC, file.path(out_main, "fig1c_proteomics_comparability"), w = 6, h = 4)

# ------------------------------------------------------------------------------
# PANEL D — Lipidomics depth (Intensity > 1e5), POS/NEG
#   solvent-first order (as for proteomics): BAMM/IPA/MeOH/standard, 10 then 5 µL, QC last
# ------------------------------------------------------------------------------
metaP <- read_tsv(meta_pos, show_col_types = FALSE) %>%
  mutate(Condition = norm_ul(Condition))
metaN <- read_tsv(meta_neg, show_col_types = FALSE) %>%
  mutate(Condition = norm_ul(Condition))

matP <- read_tsv(in_pos, show_col_types = FALSE) %>%
  mutate(.feature_id = row_number())
matN <- read_tsv(in_neg, show_col_types = FALSE) %>%
  mutate(.feature_id = row_number())

sampP <- intersect(colnames(matP), metaP$Sample_ID)
sampN <- intersect(colnames(matN), metaN$Sample_ID)

longP <- matP %>%
  select(all_of(sampP), .feature_id) %>%
  pivot_longer(cols = all_of(sampP), names_to = "Sample_ID", values_to = "Intensity") %>%
  left_join(metaP, by = "Sample_ID") %>%
  mutate(Polarity = "Positive polarity")

longN <- matN %>%
  select(all_of(sampN), .feature_id) %>%
  pivot_longer(cols = all_of(sampN), names_to = "Sample_ID", values_to = "Intensity") %>%
  left_join(metaN, by = "Sample_ID") %>%
  mutate(Polarity = "Negative polarity")

lip_all <- bind_rows(longP, longN) %>%
  mutate(
    Polarity  = factor(Polarity, levels = c("Negative polarity","Positive polarity")),
    Intensity = as.numeric(Intensity),
    Condition = norm_ul(as.character(Condition))
  )

depth_by_sample_lip <- lip_all %>%
  group_by(Polarity, Sample_ID, Condition) %>%
  summarise(n_detected = sum(is_detected(Intensity)), .groups = "drop") %>%
  filter(Condition != "BLK")

depth_sum_lip <- depth_by_sample_lip %>%
  group_by(Polarity, Condition) %>%
  summarise(Mean = mean(n_detected), SD = sd(n_detected), .groups = "drop") %>%
  mutate(
    solvent = case_when(
      Condition == "QC" ~ "QC",
      str_detect(Condition, "standard") ~ "standard",
      str_detect(Condition, "BAMM") ~ "BAMM",
      str_detect(Condition, "IPA")  ~ "IPA",
      str_detect(Condition, "MeOH") ~ "MeOH",
      TRUE ~ "Other"
    ),
    volume = case_when(
      str_detect(Condition, "10\\s*µL") ~ "10",
      str_detect(Condition, "5\\s*µL")  ~ "5",
      TRUE ~ NA_character_
    )
  )

ord_lip <- depth_sum_lip %>%
  filter(!Condition %in% c("QC", "BLK")) %>%
  distinct(Condition, solvent, volume) %>%
  mutate(
    solvent = factor(solvent, levels = c("BAMM","IPA","MeOH","standard","Other")),
    volume  = factor(volume,  levels = c("10","5"))
  ) %>%
  arrange(solvent, volume) %>%
  pull(Condition) %>%
  unique()

if ("QC" %in% depth_sum_lip$Condition) ord_lip <- c(ord_lip, "QC")


depth_sum_lip <- depth_sum_lip %>%
  mutate(CondPanel = paste0(as.character(Polarity), "__", Condition),
         CondPanel = factor(CondPanel, levels = unlist(lapply(levels(lip_all$Polarity), function(polATL){
           paste0(polATL, "__", ord_lip)
         }))))

depth_by_sample_lip <- depth_by_sample_lip %>%
  mutate(CondPanel = paste0(as.character(Polarity), "__", Condition),
         CondPanel = factor(CondPanel, levels = levels(depth_sum_lip$CondPanel)))

pD <- ggplot(depth_sum_lip, aes(x = CondPanel, y = Mean, fill = sub("^.*__", "", CondPanel))) +
  geom_col(width = 0.78, alpha = 0.95, color = "grey25", linewidth = 0.25) +
  geom_errorbar(aes(ymin = pmax(0, Mean - SD), ymax = Mean + SD),
                width = 0.22, linewidth = 0.30, color = "grey45") +
  geom_point(
    data = depth_by_sample_lip,
    aes(x = CondPanel, y = n_detected),
    position = position_jitter(width = 0.12, height = 0),
    size = 2.0, alpha = 0.30, color = "grey10",
    inherit.aes = FALSE
  ) +
  geom_text(
    aes(label = as.integer(round(Mean, 0))),
    angle = 90, vjust = 0.5, hjust = 1.55,
    size = 3.2, color = "grey10"
  ) +
  facet_wrap(~Polarity, nrow = 1, scales = "free_x") +
  scale_x_discrete(labels = function(x) sub("^.*__", "", x)) +
  scale_fill_manual(values = cond_colors, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(title = expression("Lipidomics depth (Intensity > 10\u2075)"),
       x = NULL, y = "Mean IDs") +
  theme_qomics() +
  theme(plot.margin = margin(t = 15, r = 20, b = 20, l = 20),
        legend.position = "none",
        axis.text.x = element_text(angle = 50, hjust = 1)) +
  coord_cartesian(clip = "off")

pD

save_svg_png(pD, file.path(out_main, "fig1e_lipidomics_depth"), w = 5, h = 4)

# ------------------------------------------------------------------------------
# PANEL E — Lipid class composition (macro-classes), POS/NEG
#   100% stacked bars per Condition, separated by Polarity
# ------------------------------------------------------------------------------
feat_map_P <- matP %>%
  transmute(
    .feature_id = as.integer(.feature_id),
    Category    = str_squish(as.character(Category)),
    Class       = str_squish(as.character(Class))
  ) %>%
  mutate(
    Category = if_else(is.na(Category) | Category == "", "Other", Category),
    Class    = if_else(is.na(Class)    | Class    == "", "Other", Class)
  )

feat_map_N <- matN %>%
  transmute(
    .feature_id = as.integer(.feature_id),
    Category    = str_squish(as.character(Category)),
    Class       = str_squish(as.character(Class))
  ) %>%
  mutate(
    Category = if_else(is.na(Category) | Category == "", "Other", Category),
    Class    = if_else(is.na(Class)    | Class    == "", "Other", Class)
  )

det_lip <- lip_all %>%
  mutate(
    .feature_id = as.integer(.feature_id),
    Detected    = is_detected(Intensity)
  ) %>%
  filter(Detected, !Condition %in% c("BLK")) %>%
  select(Polarity, Condition, Sample_ID, .feature_id) %>%
  mutate(
    Category = case_when(
      Polarity == "Positive polarity" ~ feat_map_P$Category[match(.feature_id, feat_map_P$.feature_id)],
      Polarity == "Negative polarity" ~ feat_map_N$Category[match(.feature_id, feat_map_N$.feature_id)],
      TRUE ~ "Other"
    ),
    Class = case_when(
      Polarity == "Positive polarity" ~ feat_map_P$Class[match(.feature_id, feat_map_P$.feature_id)],
      Polarity == "Negative polarity" ~ feat_map_N$Class[match(.feature_id, feat_map_N$.feature_id)],
      TRUE ~ "Other"
    )
  ) %>%
  mutate(
    Category = if_else(is.na(Category) | Category == "", "Other", Category),
    Class    = if_else(is.na(Class)    | Class    == "", "Other", Class)
  )

# 1) start from detected lipids with Category/Class (det_lip, before rare-class collapse)
E_in <- det_lip %>% 
  distinct(Polarity, Condition, .feature_id, Category) %>%
  mutate(
    Category = str_to_lower(str_squish(Category)),
    MacroClass = case_when(
      str_detect(Category, "glycerophospholipid") ~ "Glycerophospholipids",
      str_detect(Category, "sphingolipid")        ~ "Sphingolipids",
      str_detect(Category, "glycerolipid")        ~ "Glycerolipids",
      str_detect(Category, "fatty acyl")          ~ "Fatty acyls",
      str_detect(Category, "sterol")             ~ "Sterol lipids",
      TRUE                                       ~ "Other"
    ),
    Condition = factor(as.character(Condition), levels = ord_lip),
    MacroClass = factor(
      MacroClass,
      levels = c("Glycerophospholipids","Sphingolipids","Glycerolipids","Fatty acyls","Sterol lipids","Other")
    )
  )

# 2) composition table (proportions)
plot_df_E2 <- E_in %>%
  count(Polarity, Condition, MacroClass, name = "n_feat") %>%
  group_by(Polarity, Condition) %>%
  mutate(
    prop = n_feat / sum(n_feat),
    label = if_else(prop >= 0.08, scales::percent(prop, accuracy = 1), "") # optional labels
  ) %>%
  ungroup()

# 3) palette for macro-classes
macro_cols <- c(
  "Glycerophospholipids" = "#C7DDE8",
  "Sphingolipids"        = "#D2C6D9",
  "Glycerolipids"        = "#B88A92",
  "Fatty acyls"          = "#6F8F8B",
  "Sterol lipids"        = "#B8B8B8",
  "Other"                = "#EFEFEF"
)

pE <- ggplot(plot_df_E2, aes(x = Condition, y = prop, fill = MacroClass)) +
  geom_col(width = 0.86, color = "white", linewidth = 0.5, position = "fill") +
  facet_wrap(~Polarity, nrow = 1) +
  scale_fill_manual(values = macro_cols, name = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Lipid category composition",
    x = NULL,
    y = "Composition (%)"
  ) +
  theme_qomics() +
  theme(
    plot.margin = margin(t = 15, r = 20, b = 20, l = 20),
    axis.text.x = element_text(angle = 50, hjust = 1),
    legend.position = "right"
  )

pE
save_svg_png(pE, file.path(out_main, "figS2c_lipid_class_composition"), w = 7, h = 4)


# ------------------------------------------------------------------------------
# PANEL F — Lipidomics reproducibility (CV distributions; RIDGE), POS/NEG
#   - QC removed (features already QC-filtered upstream)
#   - CV computed on RAW intensities, only detected values (>1e5)
#   - Requires >=2 detected replicates for a lipid within condition
#   - Same solvent hue; 5 µL lighter via alpha
#   - Adds reference lines at CV=20% and 30%
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggridges)
})

# fixed order (no QC)
cond_order_noQC <- c(
  "Mag-Net BAMM 10µL","Mag-Net BAMM 5µL",
  "Mag-Net IPA 10µL","Mag-Net IPA 5µL",
  "Mag-Net MeOH 10µL","Mag-Net MeOH 5µL"
)

# group aesthetics (same hue; 5µL lighter)
fill_cols_lip <- c(
  "BAMM_10µL" = "#9BBFD4", "BAMM_5µL" = "#9BBFD4",
  "IPA_10µL"  = "#C3B2D1", "IPA_5µL"  = "#C3B2D1",
  "MeOH_10µL" = "#E6D3A3", "MeOH_5µL" = "#E6D3A3"
)
alpha_cols_lip <- c(
  "BAMM_10µL" = 0.85, "BAMM_5µL" = 0.45,
  "IPA_10µL"  = 0.85, "IPA_5µL"  = 0.45,
  "MeOH_10µL" = 0.85, "MeOH_5µL" = 0.45
)

# ---- CV table ----
cv_tbl_lip <- lip_all %>%
  mutate(
    Condition = norm_ul(as.character(Condition)),
    Polarity  = factor(Polarity, levels = c("Negative polarity", "Positive polarity")),
    Detected  = is_detected(Intensity)
  ) %>%
  filter(!Condition %in% c("QC", "BLK")) %>%                     # drop QC + blanks
  filter(Condition %in% cond_order_noQC) %>%                     # keep only method conditions
  group_by(Polarity, Condition, LipidName) %>%
  summarise(
    n_det = sum(Detected, na.rm = TRUE),
    mu    = mean(Intensity[Detected], na.rm = TRUE),
    sd    = sd(Intensity[Detected], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    CV = if_else(n_det >= 2 & is.finite(mu) & mu > 0, 100 * sd / mu, NA_real_)
  ) %>%
  filter(is.finite(CV), CV >= 0) %>%
  mutate(
    Condition = factor(Condition, levels = cond_order_noQC),
    Solvent   = parse_solvent(as.character(Condition)),
    Volume    = case_when(
      str_detect(as.character(Condition), "10\\s*µL") ~ "10µL",
      str_detect(as.character(Condition), "5\\s*µL")  ~ "5µL",
      TRUE ~ NA_character_
    ),
    Group = paste0(Solvent, "_", Volume),
    Group = factor(Group, levels = c("BAMM_10µL","BAMM_5µL","IPA_10µL","IPA_5µL","MeOH_10µL","MeOH_5µL"))
  ) %>%
  filter(!is.na(Volume), !is.na(Solvent))

# ---- plot ----
pF <- ggplot(cv_tbl_lip, aes(x = CV, y = Condition, fill = Group, alpha = Group)) +
  geom_vline(xintercept = c(20, 30), linetype = 2, linewidth = 0.35, color = "grey35") +
  ggridges::geom_density_ridges(
    scale = 1.8,
    rel_min_height = 0.01,
    color = "white",
    linewidth = 0.2
  ) +
  facet_wrap(~Polarity, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = fill_cols_lip, drop = FALSE) +
  scale_alpha_manual(values = alpha_cols_lip, guide = "none") +
  coord_cartesian(xlim = c(0, 80)) +
  labs(
    title = "Lipidomics reproducibility (CV distributions)",
    x = "CV (%)",
    y = NULL
  ) +
  theme_qomics() +
  theme(
    plot.margin = margin(t = 15, r = 20, b = 20, l = 20),
    legend.position = "none",
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    strip.text = element_text(face = "bold")
  )

pF
save_svg_png(pF, file.path(out_main, "fig1f_lipidomics_cv"), w = 5.5, h = 4)
