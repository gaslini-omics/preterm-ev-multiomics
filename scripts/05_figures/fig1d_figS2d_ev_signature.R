# ==============================================================================
# Fig. 1d / Fig. S2d — EV enrichment signature (proteins + lipids), EvMap reference-based
#   left  (protein signature)          -> Fig. 1d
#   right (conserved EV lipid markers) -> Fig. S2d
# Project:   Mag-Net EV Multi-Omics — Technical validation
# Author:    Nicole Grinovero
# Date:      Feb 2026
#
# Outputs (outputs/figures/fig1_technical_validation/ev_signature/):
#   - Fig. 1d: protein EV enrichment signature (EvMap core ↑, plasma abundant ↓)
#   - Fig. S2d: conserved EV lipid markers across conditions (z-scored heatmap)
#   - source data tables (TABLE_F_*.tsv, SOURCEDATA_F_*.tsv)
#
# Requires:
#   Proteomics:
#     - Spectronaut Report.tsv with *.raw.PG.Quantity columns
#     - EvMap2025 conserved EV protein panel file (with GeneName + MarkerType)
#   Lipidomics:
#     - cleaned tech-rep tables (POS/NEG; output of 01_techrep_lip_cleaning.R) + sample annotation
#     - Saved mapping tables:
#         MAP_reference_to_MSDIAL_POS*.tsv and MAP_reference_to_MSDIAL_NEG*.tsv
#         columns include LipidName_msdial (matched MS-DIAL name)
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(stringr)
  library(readr)
  library(svglite)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# 0) PATHS
# ------------------------------------------------------------------------------

# Proteomics
report_path <- here("data","techrep","proteomics","20260206_182115_plasma_EV_multiomics_techrep_Report.tsv")

prot_panel_path <- here("data","reference","evmap2025_conservedEV_proteins_182.txt")

# Lipidomics (processed + meta)
lip_pos_path  <- here("data","techrep","lipidomics","processed","techrep_lipids_cleaned_C18_pos.txt")
lip_neg_path  <- here("data","techrep","lipidomics","processed","techrep_lipids_cleaned_C18_neg.txt")
meta_pos_path <- here("data","techrep","lipidomics","sample_meta_C18_pos.txt")
meta_neg_path <- here("data","techrep","lipidomics","sample_meta_C18_neg.txt")

# Reference-to-MS-DIAL lipid name maps (conserved EV lipids, Rai et al. 2025)
map_pos_path <- here("data","reference","MAP_reference_to_MSDIAL_POS.tsv")
map_neg_path <- here("data","reference","MAP_reference_to_MSDIAL_NEG.tsv")

# Output
out_dir <- here("outputs","figures","fig1_technical_validation","ev_signature")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1) SETTINGS
# ------------------------------------------------------------------------------
# Two conditions compared in Fig. 1d
cond_A <- "Mag-Net BAMM 10µL"
cond_B <- "Mag-Net standard 10µL"
cmp_conds <- c(cond_A, cond_B)

# Protein: top N enriched EV-core and top N depleted plasma abundant
top_n_each <- 15

primary_ev_markers <- c(
  "PDCD6IP", "FLOT1","FLOT2","NCAM1","ADAM10"
)

plasma_abundant <- c(
  "ALB","TF","A2M","SERPINA1","ORM1","ORM2","HP","HPR","AHSG",
  "FGA","FGB","FGG",
  "C3","C4A","C4B",
  "HBA1","HBB",
  "APOA1","APOA2","APOB","APOC1","APOD"
)

# Lipid detection rule (presence/absence for the heatmap)
lip_presence_threshold <- 0  # >0 works for processed tables; raise if needed (e.g., 1e4)

# Colors
COL_ENRICH <- "#b88a92ff" 
COL_DEPL   <- "#496C82"   

# ------------------------------------------------------------------------------
# 2) HELPERS
# ------------------------------------------------------------------------------
save_svg_png <- function(p, prefix, w, h, dpi = 320) {
  ggsave(paste0(prefix, ".svg"), p, width = w, height = h, device = svglite::svglite)
  ggsave(paste0(prefix, ".png"), p, width = w, height = h, dpi = dpi, bg = "white")
}

theme_qomics <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color = "grey10"),
      axis.title = element_text(color = "grey10"),
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey92", color = NA),
      plot.title = element_text(size = base_size + 2, face = "bold"),
      plot.margin = margin(6, 10, 6, 6)
    )
}

is_det <- function(x) is.finite(x) & !is.na(x) & x > 0
strip_prefix <- function(x) str_replace(x, "^\\[[0-9]+\\]\\s*", "")
parse_rep <- function(x) str_extract(x, "C\\d+")

norm_ul <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("μL", "µL") %>%
    str_replace_all("(?i)uL", "µL") %>%
    str_replace_all("\\s+", " ") %>%
    str_squish()
}

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

parse_condition_from_sample <- function(sample_col) {
  x <- strip_prefix(sample_col)
  x <- str_replace(x, "\\.raw\\.PG\\.Quantity$", "")
  x <- str_replace(x, "_C\\d+_.*$", "")
  recode_condition(x)
}

# robust extraction of the "sp" gene field (Spectronaut multi-db)
extract_id <- function(db_str, id_str, target_db = "sp") {
  if (is.na(db_str) || is.na(id_str)) return(NA_character_)
  db_values <- str_split(db_str, ";", simplify = TRUE) |> str_trim()
  id_values <- str_split(id_str, ";", simplify = TRUE) |> str_trim()
  if (length(id_values) == 1) return(id_values[1])
  idx <- which(db_values == target_db)
  if (length(idx) > 0) return(id_values[idx[1]])
  id_values[1]
}

best_pg_per_gene <- function(df_long) {
  df_long %>%
    group_by(GeneName, PG.ProteinGroups) %>%
    summarise(
      n_obs = sum(is_det(Intensity)),
      mu    = mean(Intensity, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(GeneName, desc(n_obs), desc(mu)) %>%
    group_by(GeneName) %>%
    slice(1) %>%
    ungroup() %>%
    select(GeneName, PG.ProteinGroups)
}

calc_floor_global <- function(df_long) {
  f <- df_long %>%
    filter(Condition == "Neat", is_det(Intensity)) %>%
    summarise(f = unname(quantile(Intensity, 0.01, na.rm = TRUE))) %>%
    pull(f)
  if (!is.finite(f) || is.na(f)) {
    f <- df_long %>%
      filter(is_det(Intensity)) %>%
      summarise(f = unname(quantile(Intensity, 0.001, na.rm = TRUE))) %>%
      pull(f)
  }
  f
}

gene_fc_vs_neat <- function(df_long, keep_marker_type = FALSE) {
  floor_global <- calc_floor_global(df_long)
  
  if (keep_marker_type) {
    cond_mean <- df_long %>%
      filter(is_det(Intensity)) %>%
      group_by(GeneName, MarkerType, Condition) %>%
      summarise(MeanIntensity = mean(Intensity, na.rm = TRUE), .groups = "drop")
  } else {
    cond_mean <- df_long %>%
      filter(is_det(Intensity)) %>%
      group_by(GeneName, Condition) %>%
      summarise(MeanIntensity = mean(Intensity, na.rm = TRUE), .groups = "drop")
  }
  
  neat_ref <- cond_mean %>%
    filter(Condition == "Neat") %>%
    select(GeneName, neat_mean = MeanIntensity)
  
  cond_mean %>%
    filter(Condition != "Neat") %>%
    left_join(neat_ref, by = "GeneName") %>%
    mutate(
      neat_use = if_else(is.finite(neat_mean) & !is.na(neat_mean), neat_mean, floor_global),
      log2FC   = log2((MeanIntensity + 1) / (neat_use + 1))
    )
}

# ------------------------------------------------------------------------------
# 3) PROTEIN SIGNATURE (Fig. 1d)
# ------------------------------------------------------------------------------
rep <- read_delim(report_path, delim = "\t", show_col_types = FALSE)

sample_cols <- grep("\\.raw\\.PG\\.Quantity$", colnames(rep), value = TRUE)
if (length(sample_cols) == 0) stop("No '.raw.PG.Quantity' columns detected in Report.tsv")

prot_clean <- rep %>%
  mutate(
    GeneName = if (all(c("PG.Database", "PG.Genes") %in% names(.))) {
      map2_chr(.data$PG.Database, .data$PG.Genes, extract_id)
    } else NA_character_
  )

# EvMap protein panel (needs GeneName + MarkerType)
tmp <- read_tsv(prot_panel_path, show_col_types = FALSE)
if (!("GeneName" %in% names(tmp))) stop("EvMap protein panel must have a 'GeneName' column.")
if (!("MarkerType" %in% names(tmp))) stop("EvMap protein panel must have a 'MarkerType' column (MISEV/In vitro/Other).")

prot_ev_panel <- tmp %>%
  filter(!is.na(GeneName), GeneName != "") %>%
  mutate(
    MarkerType = case_when(
      MarkerType == "MISEV"    ~ "MISEV",
      MarkerType == "In vitro" ~ "In vitro",
      TRUE                      ~ "Other"
    )
  ) %>%
  distinct(GeneName, MarkerType)

expr_long <- prot_clean %>%
  select(PG.ProteinGroups, GeneName, all_of(sample_cols)) %>%
  pivot_longer(all_of(sample_cols), names_to = "sample_col", values_to = "Intensity") %>%
  mutate(
    Condition = parse_condition_from_sample(sample_col),
    Replicate = parse_rep(sample_col)
  )

# best PG per gene (EvMap genes only)
expr_evmap0 <- expr_long %>%
  filter(!is.na(GeneName), GeneName != "") %>%
  inner_join(prot_ev_panel, by = "GeneName")

best_pg_evmap <- best_pg_per_gene(expr_evmap0)

expr_evmap <- expr_evmap0 %>%
  inner_join(best_pg_evmap, by = c("GeneName", "PG.ProteinGroups"))

# detection filter: at least 3/4 tech reps in >=1 condition
keep_evmap <- expr_evmap %>%
  group_by(GeneName, Condition) %>%
  summarise(n_det = sum(is_det(Intensity)), .groups = "drop") %>%
  group_by(GeneName) %>%
  summarise(max_det = max(n_det), .groups = "drop") %>%
  filter(max_det >= 3) %>%
  pull(GeneName)

expr_evmap <- expr_evmap %>% filter(GeneName %in% keep_evmap)

# ALL genes table (needed for plasma depletion list)
expr_all0 <- expr_long %>% filter(!is.na(GeneName), GeneName != "")
best_pg_all <- best_pg_per_gene(expr_all0)
expr_all <- expr_all0 %>% inner_join(best_pg_all, by = c("GeneName", "PG.ProteinGroups"))

keep_all <- expr_all %>%
  group_by(GeneName, Condition) %>%
  summarise(n_det = sum(is_det(Intensity)), .groups = "drop") %>%
  group_by(GeneName) %>%
  summarise(max_det = max(n_det), .groups = "drop") %>%
  filter(max_det >= 3) %>%
  pull(GeneName)

expr_all <- expr_all %>% filter(GeneName %in% keep_all)

# gene-level FC vs Neat
fc_gene_evmap <- gene_fc_vs_neat(expr_evmap, keep_marker_type = TRUE)
fc_gene_all   <- gene_fc_vs_neat(expr_all,   keep_marker_type = FALSE)

# subset to the two displayed conditions
fc_gene_evmap <- fc_gene_evmap %>% filter(Condition %in% cmp_conds)
fc_gene_all   <- fc_gene_all   %>% filter(Condition %in% cmp_conds)

# ------------------------------------------------------------------------------
# Pick marker sets to display (ROBUST across BOTH conditions) + clean ordering
#   - EV core (EvMap; MISEV/In vitro): consistently enriched in BOTH conditions
#   - Plasma abundant: consistently depleted in BOTH conditions
#   - Ordering: enriched block (high -> low), then depleted block (low -> high)
# ------------------------------------------------------------------------------

# helper: wide log2FC (one row per gene, one column per condition)
fc_wide <- function(df_fc, cmp_conds) {
  df_fc %>%
    filter(Condition %in% cmp_conds) %>%
    select(GeneName, Condition, log2FC) %>%
    distinct() %>%
    pivot_wider(names_from = Condition, values_from = log2FC)
}

# A) EV core markers (EvMap): rank by worst-case enrichment across the two methods
ev_sel <- fc_gene_evmap %>%
  filter(MarkerType %in% c("MISEV", "In vitro")) %>%
  fc_wide(cmp_conds) %>%
  mutate(
    score_consistent = pmin(.data[[cond_A]], .data[[cond_B]], na.rm = TRUE),  # worst-case enrichment
    score_mean       = rowMeans(across(all_of(c(cond_A, cond_B))), na.rm = TRUE)
  ) %>%
  arrange(desc(score_consistent), desc(score_mean)) %>%
  slice_head(n = top_n_each)


# primary markers first, then the top-ranked EvMap genes; keep top_n_each
core_genes <- unique(c(primary_ev_markers, ev_sel$GeneName))
core_genes <- core_genes[seq_len(min(top_n_each, length(core_genes)))]

# B) Plasma abundant proteins: rank by worst-case depletion across the two methods
pl_sel <- fc_gene_all %>%
  filter(GeneName %in% plasma_abundant) %>%
  fc_wide(cmp_conds) %>%
  mutate(
    score_consistent = pmax(.data[[cond_A]], .data[[cond_B]], na.rm = TRUE),  # worst-case depletion
    score_mean       = rowMeans(across(all_of(c(cond_A, cond_B))), na.rm = TRUE)
  ) %>%
  arrange(score_consistent, score_mean) %>%  # more negative first
  slice_head(n = top_n_each)

plasma_genes <- pl_sel$GeneName

# C) Assemble plot table (two blocks)
sig_prot <- bind_rows(
  fc_gene_evmap %>%
    filter(MarkerType %in% c("MISEV", "In vitro"), GeneName %in% core_genes) %>%
    mutate(Group = "EV core markers (enriched)"),
  fc_gene_all %>%
    filter(GeneName %in% plasma_genes) %>%
    mutate(Group = "Abundant plasma proteins (depleted)")
) %>%
  mutate(
    Condition = factor(Condition, levels = cmp_conds),
    Group = factor(Group, levels = c("EV core markers (enriched)",
                                     "Abundant plasma proteins (depleted)"))
  )

# ------------------------------------------------------------------------------
# D) ORDER 
# ------------------------------------------------------------------------------

ord_enriched <- sig_prot %>%
  filter(Group == "EV core markers (enriched)") %>%
  group_by(GeneName) %>%
  summarise(
    score = min(log2FC, na.rm = TRUE),  # robust enrichment across methods
    mu    = mean(log2FC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(score, mu) %>%               # ascending: low first, high last (ends up on TOP)
  pull(GeneName)

# Depleted: we want near 0 at the TOP of the depleted block, and most negative at the BOTTOM.
# Since the BOTTOM corresponds to the FIRST factor level, we must order depleted as:
# most negative first -> near 0 last.
ord_depleted <- sig_prot %>%
  filter(Group == "Abundant plasma proteins (depleted)") %>%
  group_by(GeneName) %>%
  summarise(
    score = max(log2FC, na.rm = TRUE),  # robust depletion across methods (negative)
    mu    = mean(log2FC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(score, mu) %>%               # ascending: most negative first (BOTTOM), near 0 last (TOP)
  pull(GeneName)

# Combine BOTTOM -> TOP:
#   depleted block (bottom->top), then enriched block (bottom->top)
levels_bottom_to_top <- c(ord_depleted, ord_enriched)

sig_prot <- sig_prot %>%
  mutate(GeneName_ord = factor(GeneName, levels = levels_bottom_to_top))



p_prot <- ggplot(sig_prot, aes(x = GeneName_ord, y = log2FC, fill = Group)) +
  geom_hline(yintercept = 0, linewidth = 0.8, color = "grey35") +
  geom_col(width = 0.82, color = "grey25", linewidth = 0.20) +
  facet_grid(
    rows = vars(Condition),
    switch = "y"              # strips on the opposite side of y axis
  ) +
  scale_fill_manual(values = c(
    "EV core markers (enriched)"          = scales::alpha(COL_ENRICH, 0.9),
    "Abundant plasma proteins (depleted)" = scales::alpha(COL_DEPL, 0.7)
  )) +
  labs(
    x = NULL,
    y = expression(log[2]~"(EV / Neat)")
  ) +
  theme_qomics() +
  theme(
    # put facet strips "outside" and on the right
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text.y.right = element_text(face = "bold"),
    
    # rotate gene labels
    axis.text.x = element_text(angle = 50, hjust = 1),
    
    # margins
    plot.margin = margin(t = 15, r = 20, b = 20, l = 20),
    
    legend.position = "top"
  ) +
  coord_cartesian(clip = "off")



p_prot

save_svg_png(
  p_prot,
  file.path(out_dir, "fig1d_protein_ev_signature"),
  w = 6.5, h = 5
)

# ------------------------------------------------------------------------------
# 4) LIPIDS (Fig. S2d) — EV lipid marker coherence across ALL conditions
#     - Uses ONLY saved MAPs (POS/NEG) to define the marker set
#     - Computes mean log2(Intensity+1) on detected replicates
#     - Row-wise Z-score across conditions (so conditions are directly comparable)
#     - Main heatmap + a thin right-side strip indicating POS/NEG (no facet distortion)
#     - Optional: white dots mark incomplete detection (n_det < n_rep)
# ------------------------------------------------------------------------------

# -----------------------
# 4.0) Load saved MAPs (marker definition)
# -----------------------
map_pos <- read_tsv(map_pos_path, show_col_types = FALSE)
map_neg <- read_tsv(map_neg_path, show_col_types = FALSE)

markers_pos_found <- map_pos %>%
  filter(!is.na(LipidName_msdial), LipidName_msdial != "") %>%
  pull(LipidName_msdial) %>%
  unique()

markers_neg_found <- map_neg %>%
  filter(!is.na(LipidName_msdial), LipidName_msdial != "") %>%
  pull(LipidName_msdial) %>%
  unique()

message("Lipid markers from MAPs: POS=", length(markers_pos_found),
        " | NEG=", length(markers_neg_found))

# -----------------------
# 4.1) Load MS-DIAL processed tables + metadata
# -----------------------
lip_pos <- read_tsv(lip_pos_path, show_col_types = FALSE) %>%
  mutate(
    `Alignment ID`    = as.character(`Alignment ID`),
    `Average Mz`      = as.character(`Average Mz`),
    `Average Rt(min)` = as.character(`Average Rt(min)`)
  )

lip_neg <- read_tsv(lip_neg_path, show_col_types = FALSE) %>%
  mutate(
    `Alignment ID`    = as.character(`Alignment ID`),
    `Average Mz`      = as.character(`Average Mz`),
    `Average Rt(min)` = as.character(`Average Rt(min)`)
  )

meta_pos <- read_tsv(meta_pos_path, show_col_types = FALSE)
meta_neg <- read_tsv(meta_neg_path, show_col_types = FALSE)

meta_pos <- meta_pos %>% mutate(Condition = norm_ul(Condition))
meta_neg <- meta_neg %>% mutate(Condition = norm_ul(Condition))

# -----------------------
# 4.2) Condition order (exclude QC/BLK by default)
# -----------------------
COND_ORDER <- c(
  "Mag-Net BAMM 10µL","Mag-Net BAMM 5µL",
  "Mag-Net IPA 10µL","Mag-Net IPA 5µL",
  "Mag-Net MeOH 10µL","Mag-Net MeOH 5µL",
  "Mag-Net standard 10µL","Mag-Net standard 5µL"
)

# -----------------------
# 4.3) Build long tables and merge POS+NEG
# -----------------------
lip_pos_long <- lip_pos %>%
  pivot_longer(
    cols = all_of(meta_pos$Sample_ID),
    names_to = "Sample_ID",
    values_to = "Intensity"
  ) %>%
  left_join(meta_pos, by = "Sample_ID") %>%
  mutate(
    Condition = norm_ul(Condition),
    Polarity  = "POS"
  )

lip_neg_long <- lip_neg %>%
  pivot_longer(
    cols = all_of(meta_neg$Sample_ID),
    names_to = "Sample_ID",
    values_to = "Intensity"
  ) %>%
  left_join(meta_neg, by = "Sample_ID") %>%
  mutate(
    Condition = norm_ul(Condition),
    Polarity  = "NEG"
  )

lip_long <- bind_rows(lip_pos_long, lip_neg_long) %>%
  filter(!is.na(Condition), Condition %in% COND_ORDER) %>%
  mutate(
    Condition = factor(Condition, levels = COND_ORDER),
    Polarity  = factor(Polarity, levels = c("POS","NEG")),
    Present   = is.finite(Intensity) & !is.na(Intensity) & Intensity > lip_presence_threshold,
    log2I     = log2(Intensity + 1)
  )

# -----------------------
# 4.4) Restrict to mapped EV lipid markers (from saved MAPs only)
# -----------------------
lip_markers <- lip_long %>%
  filter(
    (Polarity == "POS" & LipidName %in% markers_pos_found) |
      (Polarity == "NEG" & LipidName %in% markers_neg_found)
  )

if (nrow(lip_markers) == 0) {
  stop("No lipid markers matched. Check that LipidName strings in MS-DIAL tables match LipidName_msdial in MAPs.")
}

# -----------------------
# 4.5) Summarise: mean log2I on detected replicates + detection completeness
# -----------------------
lip_sum <- lip_markers %>%
  group_by(Polarity, LipidName, Condition) %>%
  summarise(
    n_rep = n(),
    n_det = sum(Present, na.rm = TRUE),
    det_pct = 100 * n_det / n_rep,
    mean_log2I = if_else(
      n_det > 0,
      mean(log2I[Present], na.rm = TRUE),
      NA_real_
    ),
    .groups = "drop"
  )

# -----------------------
# 4.6) Row-wise Z-score across conditions (per lipid)
#      - uses mean_log2I (NA when never detected)
#      - if a lipid is constant or has too few values, z=0
# -----------------------
lip_z <- lip_sum %>%
  group_by(Polarity, LipidName) %>%
  mutate(
    mu = mean(mean_log2I, na.rm = TRUE),
    sdv = sd(mean_log2I, na.rm = TRUE),
    z  = (mean_log2I - mu) / sdv
  ) %>%
  ungroup() %>%
  mutate(
    z = if_else(is.na(z) | is.infinite(z), 0, z)
  )

# Order lipids: POS first, then NEG; within polarity by average abundance
lip_levels <- lip_z %>%
  group_by(Polarity, LipidName) %>%
  summarise(mu = mean(mean_log2I, na.rm = TRUE), .groups = "drop") %>%
  arrange(Polarity, desc(mu)) %>%
  pull(LipidName)

lip_z <- lip_z %>%
  mutate(LipidName_ord = factor(LipidName, levels = rev(lip_levels)))

# Robust symmetric z limits
z_max <- max(abs(lip_z$z), na.rm = TRUE)
z_max <- max(1, z_max)

# -----------------------
# 4.7) Plot: main heatmap + thin POS/NEG strip on the right
# -----------------------
pol_cols <- c(POS = "#C9DCE6", NEG = "#E8D6D9")

p_main <- ggplot(lip_z, aes(x = Condition, y = LipidName_ord, fill = z)) +
  geom_tile(color = "white", linewidth = 0.25) +
  # optional: incomplete detection overlay
  geom_point(
    data = lip_z %>% filter(n_det < n_rep),
    aes(x = Condition, y = LipidName_ord),
    inherit.aes = FALSE,
    shape = 21, size = 1.6, stroke = 0.25,
    fill = "white", color = "grey25"
  ) +
  scale_fill_gradient2(
    low = alpha("#2F7F7B", 0.8),
    mid = "white",
    high = "#A15D5D",
    midpoint = 0,
    limits = c(-z_max, z_max),
    oob = scales::squish
  ) +
  labs(
    x = NULL, y = NULL,
    fill = "Z-score"
  ) +
  theme_qomics() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

# --- side strip (POS/NEG) with labels ---
pol_cols <- c(POS = "#C9DCE6", NEG = "#E8D6D9")

# 1) build strip tiles
p_strip <- ggplot(lip_z, aes(x = 1, y = LipidName_ord, fill = Polarity)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_manual(values = pol_cols, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  theme_void() +
  theme(legend.position = "none")

# 2) add one label per polarity block (POS and NEG)
#    compute y positions at the mid-point of each block (in factor coordinates)
lab_df <- lip_z %>%
  distinct(LipidName_ord, Polarity) %>%
  mutate(y_num = as.numeric(LipidName_ord)) %>%
  group_by(Polarity) %>%
  summarise(y = mean(y_num), .groups = "drop") %>%
  mutate(x = 1.12, label = as.character(Polarity))

p_strip <- p_strip +
  geom_text(
    data = lab_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    angle = 90,
    size = 3.0,
    color = "grey25",
    fontface = "bold"
  ) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5.5, 14, 5.5, 0))  # extra right room for text

# 3) combine (main heatmap already in p_main)
p_lip <- p_main + p_strip + plot_layout(widths = c(20, 1.8))
p_lip

save_svg_png(
  p_lip,
  file.path(out_dir, "figS2d_lipid_ev_markers"),
  w = 4.7, h = 3.7
)

# -----------------------
# 4.8) Save source data
# -----------------------
write_tsv(lip_z, file.path(out_dir, "TABLE_F_lipids_zscore_allConditions.tsv"))


# ------------------------------------------------------------------------------
# 5) OUTPUT TABLES (source data)
# ------------------------------------------------------------------------------

# A) Gene-level EvMap table (all EvMap genes that passed filtering)
#    Wide format: one row per gene (and MarkerType), one column per condition
fc_evmap_wide <- fc_gene_evmap %>%
  mutate(Condition = norm_ul(Condition)) %>%
  select(GeneName, MarkerType, Condition, log2FC) %>%
  distinct() %>%
  pivot_wider(
    names_from = Condition,
    values_from = log2FC,
    values_fill = NA_real_
  ) %>%
  arrange(MarkerType, GeneName)

write_tsv(fc_evmap_wide, file.path(out_dir, "TABLE_F_evmap_geneLevel_log2FC_wide.tsv"))

# B) EV score per condition (core vs global)
#    - core: MISEV + In vitro
#    - global: all EvMap 182
ev_score_core <- fc_gene_evmap %>%
  filter(MarkerType %in% c("MISEV", "In vitro")) %>%
  group_by(Condition) %>%
  summarise(
    n_genes = n_distinct(GeneName),
    median_log2FC = median(log2FC, na.rm = TRUE),
    q25 = quantile(log2FC, 0.25, na.rm = TRUE),
    q75 = quantile(log2FC, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Set = "EV core (MISEV/In vitro)")

ev_score_global <- fc_gene_evmap %>%
  group_by(Condition) %>%
  summarise(
    n_genes = n_distinct(GeneName),
    median_log2FC = median(log2FC, na.rm = TRUE),
    q25 = quantile(log2FC, 0.25, na.rm = TRUE),
    q75 = quantile(log2FC, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Set = "EV global (EvMap 182)")

ev_score_tbl <- bind_rows(ev_score_core, ev_score_global) %>%
  mutate(Condition = factor(Condition, levels = c(
    "Mag-Net BAMM 10µL","Mag-Net BAMM 5µL",
    "Mag-Net IPA 10µL","Mag-Net IPA 5µL",
    "Mag-Net MeOH 10µL","Mag-Net MeOH 5µL",
    "Mag-Net standard 10µL","Mag-Net standard 5µL"
  ))) %>%
  arrange(Set, Condition)

write_tsv(ev_score_tbl, file.path(out_dir, "TABLE_F_evScore_coreVsGlobal_perCondition.tsv"))

# C) Source data for the protein signature plot (exact genes shown in panel)
write_tsv(sig_prot, file.path(out_dir, "SOURCEDATA_F_proteinSignature_plotData.tsv"))

