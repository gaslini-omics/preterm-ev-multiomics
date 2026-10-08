# ==============================================================================
# Project: Preterm EV Multi-Omics – Clinical integration (CAMERA on modules)
# Script:  03_modules_camera_all_outcomes_multiomics.R
# Author:  Nicole Grinovero
# Date:    Jan 2026
#
# Aim:
#   Test association between cross-omics NETWORK modules (M1–Mk) and clinical outcomes
#   using limma + duplicateCorrelation (block = patient_id) + camera gene-set testing,
#   separately for PROTEINS and LIPIDS (same module definitions from network).
#
# Inputs:
#   - data/proteomics/processed/proteins_norm_imputed.txt
#   - data/lipidomics/processed/lipids_norm_imputed.txt
#   - data/metadata/sample_metadata.txt
#   - outputs/tables/preterm/integration/integration_modules_membership_annotated.tsv
#
# Outputs:
#   - outputs/tables/preterm/integration/clinical_camera/modules_camera_all_outcomes_multiomics.tsv
#       (MASTER: layer × outcome × module; includes global FDR across all tests,
#        the sensitivity test with inter-feature correlation estimated from the data,
#        and group sizes as samples and infants = Supplementary Table S5)
#   - outputs/figures/fig7_clinical_associations/camera_qc/camera_modules_heatmap_multiomics.pdf (QC)
#   - outputs/tables/preterm/integration/clinical_camera/effects_<layer>_<outcome>.tsv
#   - outputs/tables/preterm/integration/clinical_camera/drivers_top2_<layer>_<outcome>.tsv
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(limma)
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# 0) Configuration (column names as in 01_integrate_prot_lip_network.R)
# ------------------------------------------------------------------------------
tp_levels <- c("T0","T1","T2","T3","T4")

PROT_ID   <- "UniProtID"
PROT_GENE <- "GeneName"
PROT_DESC <- "PG.ProteinDescriptions"

LIP_ID    <- "LipidName"
LIP_CLASS <- "Class"
LIP_CAT   <- "Category"

SAMPLE_ID <- "sample_id"
PAT_ID    <- "patient_id"
TP_VAR    <- "timepoint"
GA_VAR    <- "ga_birth"

# Clinical variables (binary)
OUTCOMES <- c(
  "any_brain_lesion",
  "hemorrhagic_lesion",
  "sepsis",
  "mechanical_ventilation"
)

# Drivers (interpretation only)
TOPN_MODULES_FOR_DRIVERS <- 2
TOPK_DRIVERS             <- 10

# Paths
prot_expr_path <- here("data", "proteomics", "processed", "proteins_norm_imputed.txt")
lip_expr_path  <- here("data", "lipidomics",  "processed", "lipids_norm_imputed.txt")
meta_path      <- here("data", "metadata", "sample_metadata.txt")

membership_path <- here("outputs", "tables", "preterm", "integration",
                        "integration_modules_membership_annotated.tsv")

out_dir <- here("outputs", "tables", "preterm", "integration", "clinical_camera")
fig_dir <- here("outputs", "figures", "fig7_clinical_associations", "camera_qc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

out_master <- file.path(out_dir, "modules_camera_all_outcomes_multiomics.tsv")
out_heat   <- file.path(fig_dir, "camera_modules_heatmap_multiomics.pdf")

# ------------------------------------------------------------------------------
# 1) Helpers
# ------------------------------------------------------------------------------
get_sample_cols <- function(df, candidates) candidates[candidates %in% names(df)]

# robust 0/1 coercion (handles 0/1 numeric, "0"/"1", "yes"/"no", "true"/"false")
as_binary01 <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x)) return(ifelse(is.na(x), NA_integer_, as.integer(x)))
  x2 <- tolower(trimws(as.character(x)))
  out <- dplyr::case_when(
    x2 %in% c("1","yes","y","true","t","lesion","positive","pos") ~ 1L,
    x2 %in% c("0","no","n","false","f","none","negative","neg") ~ 0L,
    TRUE ~ NA_integer_
  )
  out
}

# Build feature x sample matrix from long table with ID + sample columns
make_expr_mat <- function(df, id_col, samples) {
  stopifnot(id_col %in% names(df))
  stopifnot(all(samples %in% names(df)))
  m <- df %>%
    select(all_of(id_col), all_of(samples)) %>%
    distinct(.data[[id_col]], .keep_all = TRUE) %>%
    column_to_rownames(id_col) %>%
    as.matrix()
  storage.mode(m) <- "double"
  m
}

# Build camera index list (indices into expr rownames) from membership table
build_module_index <- function(expr_ids, membership_df, feature_col = "Feature", module_col = "ModuleLabel") {
  modules <- sort(unique(membership_df[[module_col]]))
  idx <- lapply(modules, function(m) which(expr_ids %in% membership_df[[feature_col]][membership_df[[module_col]] == m]))
  names(idx) <- modules
  idx <- idx[lengths(idx) > 0]
  idx
}

# ------------------------------------------------------------------------------
# 2) Load metadata + membership
# ------------------------------------------------------------------------------
stopifnot(file.exists(prot_expr_path), file.exists(lip_expr_path),
          file.exists(meta_path), file.exists(membership_path))

meta <- read_delim(meta_path, delim = "\t", show_col_types = FALSE) %>%
  mutate(
    !!SAMPLE_ID := as.character(.data[[SAMPLE_ID]]),
    !!PAT_ID    := factor(.data[[PAT_ID]]),
    !!TP_VAR    := factor(.data[[TP_VAR]], levels = tp_levels, ordered = TRUE),
    !!GA_VAR    := as.numeric(.data[[GA_VAR]])
  )

# outcome sanity
miss_meta <- setdiff(c(SAMPLE_ID, PAT_ID, TP_VAR, GA_VAR, OUTCOMES), names(meta))
if (length(miss_meta) > 0) stop("Missing columns in metadata: ", paste(miss_meta, collapse = ", "), call. = FALSE)

# module membership from 01_integrate_prot_lip_network.R (Feature, Type, ModuleLabel + annotations)
memb <- read_delim(membership_path, delim = "\t", show_col_types = FALSE)

need_memb <- c("ModuleLabel","Type","Feature")
miss_memb <- setdiff(need_memb, names(memb))
if (length(miss_memb) > 0) stop("Missing columns in membership file: ", paste(miss_memb, collapse = ", "), call. = FALSE)

memb_prot <- memb %>% filter(Type == "Protein") %>% distinct(ModuleLabel, Feature, .keep_all = TRUE)
memb_lip  <- memb %>% filter(Type == "Lipid")   %>% distinct(ModuleLabel, Feature, .keep_all = TRUE)

# ------------------------------------------------------------------------------
# 3) Load expression tables
# ------------------------------------------------------------------------------
prot_tbl <- read_delim(prot_expr_path, delim = "\t", show_col_types = FALSE)
lip_tbl  <- read_delim(lip_expr_path,  delim = "\t", show_col_types = FALSE) %>%
  # deuterated internal standards excluded, as in the network script
  filter(!str_detect(.data[[LIP_ID]], "\\(d[0-9]+\\)"))

# match paired samples (present in both, ordered as meta)
prot_samples <- get_sample_cols(prot_tbl, meta[[SAMPLE_ID]])
lip_samples  <- get_sample_cols(lip_tbl,  meta[[SAMPLE_ID]])
paired_samples <- intersect(prot_samples, lip_samples)

if (length(paired_samples) < 10) {
  stop("Too few paired samples across proteomics+lipidomics. Check sample IDs/columns.", call. = FALSE)
}

meta_p <- meta %>%
  filter(.data[[SAMPLE_ID]] %in% paired_samples) %>%
  arrange(match(.data[[SAMPLE_ID]], paired_samples))

paired_samples <- meta_p[[SAMPLE_ID]]
message("Paired samples used: ", length(paired_samples))

# expression matrices (features x samples)
expr_prot_all <- make_expr_mat(prot_tbl, PROT_ID, paired_samples)
expr_lip_all  <- make_expr_mat(lip_tbl,  LIP_ID,  paired_samples)

# keep only features that are actually in membership (modules) for each layer
prot_ids <- intersect(rownames(expr_prot_all), memb_prot$Feature)
lip_ids  <- intersect(rownames(expr_lip_all),  memb_lip$Feature)

expr_prot <- expr_prot_all[prot_ids, , drop = FALSE]
expr_lip  <- expr_lip_all[lip_ids,  , drop = FALSE]

# module indices per layer
idx_prot <- build_module_index(rownames(expr_prot), memb_prot, feature_col = "Feature", module_col = "ModuleLabel")
idx_lip  <- build_module_index(rownames(expr_lip),  memb_lip,  feature_col = "Feature", module_col = "ModuleLabel")

if (length(idx_prot) == 0) stop("No non-empty protein module sets after filtering.", call. = FALSE)
if (length(idx_lip)  == 0) stop("No non-empty lipid module sets after filtering.", call. = FALSE)

message("Modules tested (proteins): ", paste(names(idx_prot), collapse = ", "))
message("Modules tested (lipids):   ", paste(names(idx_lip),  collapse = ", "))

# ------------------------------------------------------------------------------
# 4) Core runner: limma (duplicateCorrelation) + camera for one layer + one outcome
# ------------------------------------------------------------------------------
run_layer_outcome <- function(layer, expr_mat, idx_list, meta_df, outcome_var, annot_df = NULL) {
  
  stopifnot(all(colnames(expr_mat) == meta_df[[SAMPLE_ID]]))
  if (!outcome_var %in% names(meta_df)) stop("Outcome not in metadata: ", outcome_var, call. = FALSE)
  
  # coerce outcome to 0/1
  y01 <- as_binary01(meta_df[[outcome_var]])
  keep <- !is.na(y01)
  
  if (sum(keep) < 10) {
    warning("[", layer, " | ", outcome_var, "] Too few non-missing samples. Skipping.")
    return(list(camera = tibble(), effects = tibble(), drivers = tibble()))
  }
  
  meta_o <- meta_df[keep, , drop = FALSE]
  expr_o <- expr_mat[, meta_o[[SAMPLE_ID]], drop = FALSE]
  
  meta_o$clin <- factor(y01[keep], levels = c(0,1), labels = c("No","Yes"))
  
  tab <- table(meta_o$clin)
  message("[", layer, " | ", outcome_var, "] N(No)=", tab["No"], " | N(Yes)=", tab["Yes"])

  # group sizes (samples and infants) for Supplementary Table S5, counted on the
  # full cohort in the metadata (74 samples). The tests below use only the
  # samples with both proteomic and lipidomic data (73; one sample has no
  # lipidomic profile).
  y_all <- as_binary01(meta[[outcome_var]])
  sizes <- tibble(
    n_samples_no      = sum(y_all == 0, na.rm = TRUE),
    n_samples_yes     = sum(y_all == 1, na.rm = TRUE),
    n_infants_no      = n_distinct(meta[[PAT_ID]][which(y_all == 0)]),
    n_infants_yes     = n_distinct(meta[[PAT_ID]][which(y_all == 1)]),
    n_samples_missing = sum(is.na(y_all))
  )
  
  # design: clin + timepoint + GA
  design <- model.matrix(~ clin + timepoint + ga_birth, data = meta_o)
  
  # repeated measures via duplicateCorrelation
  corfit <- duplicateCorrelation(expr_o, design, block = meta_o[[PAT_ID]])
  fit <- lmFit(expr_o, design, block = meta_o[[PAT_ID]], correlation = corfit$consensus)
  fit <- eBayes(fit)
  
  coef_name <- "clinYes"
  if (!coef_name %in% colnames(fit$coefficients)) {
    stop("[", layer, " | ", outcome_var, "] coef ", coef_name, " not found (check encoding).", call. = FALSE)
  }
  
  # feature-level effects (rich, but per outcome per layer)
  eff <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none") %>%
    rownames_to_column("Feature") %>%
    as_tibble() %>%
    transmute(
      Layer   = layer,
      Outcome = outcome_var,
      Feature = Feature,
      beta    = logFC,
      t       = t,
      p       = P.Value,
      FDR     = adj.P.Val
    ) %>%
    arrange(FDR, desc(abs(beta)))
  
  # camera set test
  contrast <- rep(0, ncol(design)); names(contrast) <- colnames(design)
  contrast[coef_name] <- 1
  
  cam <- camera(expr_o, index = idx_list, design = design, contrast = contrast) %>%
    as.data.frame() %>%
    rownames_to_column("ModuleLabel") %>%
    as_tibble() %>%
    transmute(
      Layer     = layer,
      Outcome   = outcome_var,
      ModuleLabel,
      N         = NGenes,
      Direction = Direction,   # Up = higher in Yes; Down = lower in Yes
      PValue    = PValue
    )

  # sensitivity: inter-feature correlation estimated from the data
  # (inter.gene.cor = NA) instead of the limma default (0.01)
  cam_est <- camera(expr_o, index = idx_list, design = design, contrast = contrast,
                    inter.gene.cor = NA) %>%
    as.data.frame() %>%
    rownames_to_column("ModuleLabel") %>%
    as_tibble() %>%
    transmute(ModuleLabel,
              P_estimated_igc = PValue,
              IGC_estimated   = Correlation)

  cam <- cam %>%
    left_join(cam_est, by = "ModuleLabel") %>%
    bind_cols(sizes[rep(1, nrow(cam)), ])
  
  # drivers (interpretation only): top2 modules by PValue, top10 features by p
  top_mods <- cam %>%
    arrange(PValue) %>%
    slice_head(n = TOPN_MODULES_FOR_DRIVERS) %>%
    pull(ModuleLabel)
  
  drivers <- tibble()
  if (length(top_mods) > 0) {
    memb_layer <- if (layer == "Protein") memb_prot else memb_lip
    
    drivers <- memb_layer %>%
      filter(ModuleLabel %in% top_mods) %>%
      transmute(ModuleLabel, Feature = Feature) %>%
      inner_join(eff, by = c("Feature")) %>%
      group_by(Layer, Outcome, ModuleLabel) %>%
      arrange(p, desc(abs(beta)), .by_group = TRUE) %>%
      slice_head(n = TOPK_DRIVERS) %>%
      ungroup()
    
    # add feature annotation if provided
    if (!is.null(annot_df)) {
      drivers <- drivers %>%
        left_join(annot_df, by = "Feature")
    }
  }
  
  list(camera = cam, effects = eff, drivers = drivers)
}

# Build annotation tables from membership output (already includes gene/desc for proteins, class/cat for lipids)
prot_annot <- memb_prot %>%
  transmute(
    Feature = Feature,
    GeneName = if (PROT_GENE %in% names(memb_prot)) .data[[PROT_GENE]] else NA_character_,
    Description = if (PROT_DESC %in% names(memb_prot)) .data[[PROT_DESC]] else NA_character_
  ) %>% distinct()

lip_annot <- memb_lip %>%
  transmute(
    Feature = Feature,
    Class = if (LIP_CLASS %in% names(memb_lip)) .data[[LIP_CLASS]] else NA_character_,
    Category = if (LIP_CAT %in% names(memb_lip)) .data[[LIP_CAT]] else NA_character_
  ) %>% distinct()

# ------------------------------------------------------------------------------
# 5) Run ALL outcomes for BOTH layers
# ------------------------------------------------------------------------------
all_cam <- list()
all_eff <- list()
all_drv <- list()

for (oc in OUTCOMES) {
  
  # PROTEINS
  res_p <- run_layer_outcome(
    layer = "Protein",
    expr_mat = expr_prot,
    idx_list = idx_prot,
    meta_df  = meta_p,
    outcome_var = oc,
    annot_df = prot_annot
  )
  
  write_tsv(res_p$effects,  file.path(out_dir, paste0("effects_Protein_", oc, ".tsv")))
  write_tsv(res_p$drivers,  file.path(out_dir, paste0("drivers_top2_Protein_", oc, ".tsv")))
  
  all_cam[[paste0("Protein__", oc)]] <- res_p$camera
  all_eff[[paste0("Protein__", oc)]] <- res_p$effects
  all_drv[[paste0("Protein__", oc)]] <- res_p$drivers
  
  # LIPIDS
  res_l <- run_layer_outcome(
    layer = "Lipid",
    expr_mat = expr_lip,
    idx_list = idx_lip,
    meta_df  = meta_p,
    outcome_var = oc,
    annot_df = lip_annot
  )
  
  write_tsv(res_l$effects,  file.path(out_dir, paste0("effects_Lipid_", oc, ".tsv")))
  write_tsv(res_l$drivers,  file.path(out_dir, paste0("drivers_top2_Lipid_", oc, ".tsv")))
  
  all_cam[[paste0("Lipid__", oc)]] <- res_l$camera
  all_eff[[paste0("Lipid__", oc)]] <- res_l$effects
  all_drv[[paste0("Lipid__", oc)]] <- res_l$drivers
}

cam_tbl <- bind_rows(all_cam)

# Global FDR across ALL tests:
# layers (2) × outcomes (4) × modules (~5) = ~40 tests
cam_tbl <- cam_tbl %>%
  mutate(
    FDR_global = p.adjust(PValue, method = "BH"),
    FDR_estimated_igc = p.adjust(P_estimated_igc, method = "BH"),
    Outcome = factor(Outcome, levels = OUTCOMES),
    ModuleLabel = factor(ModuleLabel, levels = sort(unique(ModuleLabel)))
  ) %>%
  arrange(FDR_global, Layer, Outcome, ModuleLabel) %>%
  select(Layer, Outcome, ModuleLabel, N, Direction, PValue, FDR_global,
         P_estimated_igc, IGC_estimated, FDR_estimated_igc,
         n_samples_no, n_samples_yes, n_infants_no, n_infants_yes, n_samples_missing)

write_tsv(cam_tbl, out_master)

# ------------------------------------------------------------------------------
# 6) Heatmap (facet by layer): signed -log10(PValue)
# ------------------------------------------------------------------------------
heat_df <- cam_tbl %>%
  mutate(
    signed_score = ifelse(Direction == "Up", 1, -1) * (-log10(PValue))
  )

p_heat <- ggplot(heat_df, aes(x = Outcome, y = ModuleLabel, fill = signed_score)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(FDR_global < 0.05, "*", "")), size = 6) +
  facet_wrap(~ Layer, ncol = 1) +
  theme_minimal(base_size = 11) +
  labs(
    title = "Module set associations (camera) – Proteins vs Lipids",
    subtitle = "Fill = signed -log10(P) (Up = higher in Yes; Down = lower in Yes). * FDR_global < 0.05",
    x = NULL, y = NULL
  ) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

p_heat

ggsave(out_heat, p_heat, width = 9.2, height = 6.0)

cat("DONE.\nSaved:\n- ", out_master, "\n- ", out_heat, "\n", sep = "")
cat("Per-outcome effects + top2 drivers saved in:\n- ", out_dir, "\n", sep = "")
