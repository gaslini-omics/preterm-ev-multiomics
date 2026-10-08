# ==============================================================================
# Project: Preterm EV Multi-Omics – Integration (proteomics + lipidomics)
# Script:  01_integrate_prot_lip_network.R
# Author:  Nicole Grinovero
# Date:    Jan 2026 (updated Oct 2026)
#
# Aim:
#   Build a reproducible cross-omics correlation network (proteins ↔ lipids)
#   using features selected by longitudinal LMM (FDR_global < 0.05),
#   identify data-driven cross-omics modules via Louvain community detection,
#   and generate network-level outputs for downstream interpretation and figure
#   generation.
#   Assess network robustness (section 10): correlations residualized for
#   timepoint and gestational age, repeated-measures correlations, and
#   leave-one-infant-out stability of edges and modules.
#
# Inputs:
#   - data/proteomics/processed/proteins_norm_imputed.txt
#   - data/lipidomics/processed/lipids_norm_imputed.txt
#   - data/metadata/sample_metadata.txt
#   - outputs/tables/preterm/proteomics/proteins_LMM_time_results.txt
#   - outputs/tables/preterm/lipidomics/lipids_LMM_time_results.txt
#
# Outputs (outputs/tables/preterm/integration/):
#   - integration_correlation_landscape.tsv.gz
#       Full protein–lipid correlation landscape (all tested pairs)
#   - integration_correlation_landscape_stats.tsv
#       Summary statistics of cross-omics coordination 
#   - integration_edges.tsv
#       Significant protein–lipid edges (FDR < threshold)
#   - integration_edges_filtered.tsv
#       Network edges retained after module-size filtering (Cytoscape-ready)
#   - integration_nodes.tsv
#       Node attributes table (modules, degree, annotations; Cytoscape-ready)
#   - integration_modules_membership_annotated.tsv
#       Full module membership with annotations and network metrics
#   - integration_modules_summary.tsv
#       Module composition summary (proteins, lipids, total)
#   - integration_network.graphml
#       Filtered cross-omics network for visualization in Cytoscape
#   - network_robustness_edges_summary.tsv
#   - network_robustness_loio_partition.tsv
#   - network_robustness_louvain_baseline.tsv
#       Robustness of edges and modules (section 10)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(tibble)
  library(here)
  library(igraph)
})

set.seed(1)

# ------------------------------------------------------------------------------
# 0) Configuration
# ------------------------------------------------------------------------------
FDR_LMM         <- 0.05
FDR_EDGE        <- 0.05
RHO_FOR_NETWORK <- 0.35
MIN_MOD         <- 10

SAVE_GRAPHML <- TRUE

tp_levels <- c("T0","T1","T2","T3","T4")

PROT_ID   <- "UniProtID"
PROT_GENE <- "GeneName"
PROT_DESC <- "PG.ProteinDescriptions"

LIP_ID    <- "LipidName"
LIP_CLASS <- "Class"
LIP_CAT   <- "Category"

prot_expr_path <- here("data", "proteomics", "processed", "proteins_norm_imputed.txt")
lip_expr_path  <- here("data", "lipidomics",  "processed", "lipids_norm_imputed.txt")
meta_path      <- here("data", "metadata", "sample_metadata.txt")

prot_lmm_path  <- here("outputs", "tables", "preterm", "proteomics", "proteins_LMM_time_results.txt")
lip_lmm_path   <- here("outputs", "tables", "preterm", "lipidomics",  "lipids_LMM_time_results.txt")

out_dir <- here("outputs", "tables", "preterm", "integration")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1) Minimal helpers
# ------------------------------------------------------------------------------
make_mat <- function(df, id_col, sample_cols) {
  df %>%
    select(all_of(id_col), all_of(sample_cols)) %>%
    distinct(.data[[id_col]], .keep_all = TRUE) %>%
    column_to_rownames(id_col) %>%
    as.matrix() %>%
    t()
}

z_by_feature <- function(m) {
  m2 <- scale(m, center = TRUE, scale = TRUE)
  m2[is.na(m2)] <- 0
  m2
}

r_to_p <- function(r, df) {
  r <- pmin(pmax(r, -0.999999), 0.999999)
  tval <- r * sqrt(df / (1 - r^2))
  2 * pt(abs(tval), df = df, lower.tail = FALSE)
}

# ------------------------------------------------------------------------------
# 2) Load
# ------------------------------------------------------------------------------
prot <- read_delim(prot_expr_path, delim = "\t", show_col_types = FALSE)
lip  <- read_delim(lip_expr_path,  delim = "\t", show_col_types = FALSE) %>%
  filter(!str_detect(.data[[LIP_ID]], "\\(d[0-9]+\\)"))

meta <- read_delim(meta_path, delim = "\t", show_col_types = FALSE) %>%
  mutate(
    sample_id  = as.character(sample_id),
    patient_id = factor(patient_id),
    timepoint  = factor(timepoint, levels = tp_levels, ordered = TRUE)
  )

prot_lmm <- read_delim(prot_lmm_path, delim = "\t", show_col_types = FALSE)
lip_lmm  <- read_delim(lip_lmm_path,  delim = "\t", show_col_types = FALSE)

# ------------------------------------------------------------------------------
# 3) Paired samples 
# ------------------------------------------------------------------------------
paired_samples <- intersect(intersect(meta$sample_id, colnames(prot)), intersect(meta$sample_id, colnames(lip)))

meta_p <- meta %>%
  filter(sample_id %in% paired_samples) %>%
  arrange(match(sample_id, paired_samples))

paired_samples <- meta_p$sample_id
cat("Paired samples:", length(paired_samples), "\n")

# ------------------------------------------------------------------------------
# 4) LMM-selected features
# ------------------------------------------------------------------------------
prot_sig <- prot_lmm %>% filter(FDR_global < FDR_LMM) %>% pull(.data[[PROT_ID]]) %>% unique()
lip_sig  <- lip_lmm  %>% filter(FDR_global < FDR_LMM) %>% pull(.data[[LIP_ID]])  %>% unique()

prot_mat_all <- make_mat(prot, PROT_ID, paired_samples)
lip_mat_all  <- make_mat(lip,  LIP_ID,  paired_samples)

prot_keep <- intersect(colnames(prot_mat_all), prot_sig)
lip_keep  <- intersect(colnames(lip_mat_all),  lip_sig)

prot_mat <- prot_mat_all[, prot_keep, drop = FALSE]
lip_mat  <- lip_mat_all[,  lip_keep,  drop = FALSE]

cat("Temporal proteins:", ncol(prot_mat), "\n")
cat("Temporal lipids:  ", ncol(lip_mat),  "\n")

# ------------------------------------------------------------------------------
# 5) Spearman correlations + BH-FDR  
# ------------------------------------------------------------------------------
prot_z <- z_by_feature(prot_mat)
lip_z  <- z_by_feature(lip_mat)

prot_r <- apply(prot_z, 2, rank, ties.method = "average")
lip_r  <- apply(lip_z,  2, rank, ties.method = "average")

R_pl <- cor(prot_r, lip_r, method = "pearson", use = "pairwise.complete.obs")

P_pl <- matrix(
  r_to_p(as.vector(R_pl), nrow(prot_z) - 2),
  nrow = nrow(R_pl), ncol = ncol(R_pl),
  dimnames = dimnames(R_pl)
)

cor_df <- as.data.frame(as.table(R_pl)) %>%
  rename(Protein = Var1, Lipid = Var2, rho = Freq) %>%
  mutate(
    p       = as.vector(P_pl),
    FDR     = p.adjust(p, method = "BH"),
    abs_rho = abs(rho),
    sign    = ifelse(rho > 0, "pos", "neg"),
    significant = (FDR < FDR_EDGE)
  )

# Save full correlation landscape
write_tsv(cor_df, file.path(out_dir, "integration_correlation_landscape.tsv.gz"))

edges_sig_FDR <- cor_df %>% filter(significant)
prot_involved <- unique(edges_sig_FDR$Protein)
lip_involved  <- unique(edges_sig_FDR$Lipid)

stats_tbl <- tibble(
  n_samples_paired   = nrow(prot_z),
  n_prot_temporal    = ncol(prot_z),
  n_lip_temporal     = ncol(lip_z),
  n_tests_total      = ncol(prot_z) * ncol(lip_z),
  n_edges_FDR        = nrow(edges_sig_FDR),
  n_edges_FDR_pos    = sum(edges_sig_FDR$rho > 0, na.rm = TRUE),
  n_edges_FDR_neg    = sum(edges_sig_FDR$rho < 0, na.rm = TRUE),
  median_abs_rho_FDR = if (nrow(edges_sig_FDR) > 0) median(edges_sig_FDR$abs_rho, na.rm = TRUE) else NA_real_,
  min_rho_FDR        = if (nrow(edges_sig_FDR) > 0) min(edges_sig_FDR$rho, na.rm = TRUE) else NA_real_,
  max_rho_FDR        = if (nrow(edges_sig_FDR) > 0) max(edges_sig_FDR$rho, na.rm = TRUE) else NA_real_,
  n_prot_involved    = length(prot_involved),
  n_lip_involved     = length(lip_involved),
  pct_prot_involved  = 100 * length(prot_involved) / ncol(prot_z),
  pct_lip_involved   = 100 * length(lip_involved)  / ncol(lip_z)
)

write_tsv(stats_tbl, file.path(out_dir, "integration_correlation_landscape_stats.tsv"))

# ------------------------------------------------------------------------------
# 6) Network edges used for module detection (FDR + |rho|)
# ------------------------------------------------------------------------------
edges_sig <- cor_df %>%
  filter(FDR < FDR_EDGE, abs_rho >= RHO_FOR_NETWORK) %>%
  mutate(
    Protein = as.character(Protein),
    Lipid   = as.character(Lipid),
    node_min = pmin(Protein, Lipid),
    node_max = pmax(Protein, Lipid),
    weight   = abs_rho
  ) %>%
  arrange(node_min, node_max, desc(weight), FDR) %>%
  select(Protein, Lipid, rho, abs_rho, sign, weight)

cat("Edges kept for network (FDR + |rho|):", nrow(edges_sig), "\n")

write_tsv(edges_sig, file.path(out_dir, "integration_edges.tsv"))

# ------------------------------------------------------------------------------
# 7) Graph + Louvain modules 
# ------------------------------------------------------------------------------
node_tbl <- tibble(name = c(colnames(prot_z), colnames(lip_z))) %>%
  mutate(Type = ifelse(name %in% colnames(prot_z), "Protein", "Lipid")) %>%
  arrange(Type, name) 

g <- graph_from_data_frame(
  d = edges_sig %>% transmute(from = Protein, to = Lipid, rho, abs_rho, sign, weight),
  vertices = node_tbl,
  directed = FALSE
)

# Louvain
cl <- cluster_louvain(g, weights = E(g)$weight)
V(g)$Module_old <- as.integer(membership(cl))

mod_sum <- tibble(Feature = V(g)$name, Type = V(g)$Type, Module_old = V(g)$Module_old) %>%
  count(Module_old, Type, name = "n") %>%
  pivot_wider(names_from = Type, values_from = n, values_fill = 0) %>%
  mutate(Total = Protein + Lipid) %>%
  arrange(desc(Total), Module_old)

keep_old <- mod_sum %>% filter(Total >= MIN_MOD) %>% pull(Module_old)

# deterministic relabel: by size (desc), ties by Module_old
mod_map <- mod_sum %>%
  filter(Module_old %in% keep_old) %>%
  arrange(desc(Total), Module_old) %>%
  mutate(Module = row_number(),
         ModuleLabel = paste0("M", Module)) %>%
  select(Module_old, Module, ModuleLabel, Protein, Lipid, Total)

write_tsv(
  mod_map %>% select(Module, ModuleLabel, Protein, Lipid, Total) %>% arrange(Module),
  file.path(out_dir, "integration_modules_summary.tsv")
)

g_f <- induced_subgraph(g, vids = V(g)[Module_old %in% keep_old])
V(g_f)$Module <- mod_map$Module[match(V(g_f)$Module_old, mod_map$Module_old)]
V(g_f)$ModuleLabel <- paste0("M", V(g_f)$Module)

# Export edges AFTER module-size filter (Cytoscape-ready, deterministic sort)
edges_f <- as_data_frame(g_f, what = "edges") %>%
  transmute(
    Protein = from,
    Lipid   = to,
    rho     = rho,
    abs_rho = abs_rho,
    sign    = sign,
    weight  = abs_rho,
    node_min = pmin(from, to),
    node_max = pmax(from, to)
  ) %>%
  arrange(node_min, node_max, desc(weight)) %>%
  select(-node_min, -node_max)

write_tsv(edges_f, file.path(out_dir, "integration_edges_filtered.tsv"))
cat("Edges kept (post-module filter):", nrow(edges_f), "\n")

# ------------------------------------------------------------------------------
# 8) Membership + annotations + node metrics + nodes.tsv 
# ------------------------------------------------------------------------------
prot_annot <- prot %>%
  select(any_of(c(PROT_ID, PROT_GENE, PROT_DESC))) %>%
  distinct() %>%
  rename(Feature = all_of(PROT_ID)) %>%
  mutate(Type = "Protein")

lip_annot <- lip %>%
  select(any_of(c(LIP_ID, LIP_CLASS, LIP_CAT))) %>%
  distinct() %>%
  rename(Feature = all_of(LIP_ID)) %>%
  mutate(Type = "Lipid")

annot_tbl <- bind_rows(prot_annot, lip_annot)

deg_tbl <- tibble(Feature = V(g_f)$name, Degree = degree(g_f))

rho_node <- bind_rows(
  edges_f %>% transmute(Feature = Protein, rho, abs_rho),
  edges_f %>% transmute(Feature = Lipid,   rho, abs_rho)
) %>%
  group_by(Feature) %>%
  summarise(
    MeanAbsRho = mean(abs_rho, na.rm = TRUE),
    MeanRho    = mean(rho, na.rm = TRUE),
    .groups = "drop"
  )

membership_tbl <- tibble(
  Feature = V(g_f)$name,
  Type    = V(g_f)$Type,
  Module  = V(g_f)$Module,
  ModuleLabel = V(g_f)$ModuleLabel
) %>%
  left_join(deg_tbl,   by = "Feature") %>%
  left_join(rho_node,  by = "Feature") %>%
  left_join(annot_tbl, by = c("Feature","Type")) %>%
  arrange(Module, desc(Degree), Type, Feature)

write_tsv(membership_tbl, file.path(out_dir, "integration_modules_membership_annotated.tsv"))

nodes_out <- membership_tbl %>%
  transmute(
    id = Feature,
    type = Type,
    module = ModuleLabel,
    module_n = Module,
    degree = Degree,
    mean_abs_rho = MeanAbsRho,
    mean_rho = MeanRho,
    gene = if (PROT_GENE %in% names(membership_tbl)) .data[[PROT_GENE]] else NA_character_,
    description = if (PROT_DESC %in% names(membership_tbl)) .data[[PROT_DESC]] else NA_character_,
    lipid_class = if (LIP_CLASS %in% names(membership_tbl)) .data[[LIP_CLASS]] else NA_character_,
    lipid_category = if (LIP_CAT %in% names(membership_tbl)) .data[[LIP_CAT]] else NA_character_
  )

write_tsv(nodes_out, file.path(out_dir, "integration_nodes.tsv"))

# ------------------------------------------------------------------------------
# 9) Optional GraphML
# ------------------------------------------------------------------------------
if (isTRUE(SAVE_GRAPHML)) {
  write_graph(g_f, file.path(out_dir, "integration_network.graphml"), format = "graphml")
}

cat("Done.\nSaved in: ", out_dir, "\n", sep = "")

# ==============================================================================
# 10) Network robustness
#
#   Reference = network edges (edges_sig) and Louvain modules built above.
#   (a) correlations on ranks residualized for timepoint + ga_birth
#       (partial correlation; df = n - 2 - number of covariate columns)
#   (b) repeated-measures correlation: ranks centered within infant
#       (removes all between-infant differences, incl. ga_birth; df = n - k - 1)
#   (c) leave-one-infant-out (LOIO): edges and Louvain partition rebuilt
#       without each infant in turn
#   Louvain is stochastic: partition agreement (adjusted Rand index, ARI) is
#   compared with the run-to-run variability on the full data.
#
#   Outputs (outputs/tables/preterm/integration/):
#     - network_robustness_edges_summary.tsv
#     - network_robustness_loio_partition.tsv
#     - network_robustness_louvain_baseline.tsv
# ==============================================================================
RUN_ROBUSTNESS <- TRUE
N_RUNS_FULL    <- 50   # Louvain runs on full data (baseline)
N_RUNS_LOIO    <- 20   # Louvain runs per reduced dataset

if (isTRUE(RUN_ROBUSTNESS)) {

  rob_dir <- out_dir

  infant <- droplevels(meta_p$patient_id)
  tp     <- factor(as.character(meta_p$timepoint), levels = tp_levels)
  ga     <- meta_p$ga_birth

  # correlations from rank matrices (as in section 5), three modes
  cor_from_ranks <- function(Pr, Lr, infant, tp, ga, mode = "pooled") {
    n <- nrow(Pr)
    if (mode == "within") {
      cen <- function(m) m - apply(m, 2, function(x) ave(x, infant))
      Pr <- cen(Pr); Lr <- cen(Lr); df <- n - nlevels(droplevels(infant)) - 1
    } else if (mode == "resid") {
      X <- model.matrix(~ droplevels(tp) + ga)
      res <- function(m) m - X %*% qr.solve(X, m)
      Pr <- res(Pr); Lr <- res(Lr); df <- n - 2 - (ncol(X) - 1)
    } else {
      df <- n - 2
    }
    R <- cor(Pr, Lr)
    p <- r_to_p(as.vector(R), df)
    tibble(Protein = rep(rownames(R), ncol(R)), Lipid = rep(colnames(R), each = nrow(R)),
           rho = as.vector(R), p = p, FDR = p.adjust(p, method = "BH"))
  }
  ranks_of <- function(m) apply(z_by_feature(m), 2, rank, ties.method = "average")

  # edge-level comparison with the reference network edges
  compare_edges <- function(cd, label) {
    b <- inner_join(edges_sig %>% select(Protein, Lipid, rho_ref = rho), cd, by = c("Protein", "Lipid"))
    same <- sign(b$rho) == sign(b$rho_ref)
    tibble(analysis = label, n_edges = nrow(b),
           pct_same_sign         = 100 * mean(same),
           pct_FDR_same_sign     = 100 * mean(b$FDR < FDR_EDGE & same),
           pct_nominal_same_sign = 100 * mean(b$p < 0.05 & same),
           median_abs_rho_ref    = median(abs(b$rho_ref)),
           median_abs_rho        = median(abs(b$rho)))
  }

  # (a) + (b)
  edge_tbl <- bind_rows(
    compare_edges(cor_from_ranks(prot_r, lip_r, infant, tp, ga, "resid"),  "residualized_timepoint_GA"),
    compare_edges(cor_from_ranks(prot_r, lip_r, infant, tp, ga, "within"), "within_infant_rmcorr")
  )

  # partition helpers (reference = modules of this run)
  mod_ref <- setNames(V(g_f)$Module, V(g_f)$name)
  ari <- function(a, b) {
    tab <- table(a, b); c2 <- function(x) x * (x - 1) / 2
    s <- sum(c2(tab)); sa <- sum(c2(rowSums(tab))); sb <- sum(c2(colSums(tab)))
    e <- sa * sb / c2(sum(tab)); (s - e) / (0.5 * (sa + sb) - e)
  }
  build_graph <- function(e) graph_from_data_frame(
    e %>% transmute(from = Protein, to = Lipid, weight = abs(rho)),
    vertices = node_tbl %>% select(name), directed = FALSE)
  louvain_run <- function(g, seed) {
    set.seed(seed)
    setNames(as.integer(membership(cluster_louvain(g, weights = E(g)$weight))), V(g)$name)
  }
  ari_ref  <- function(m) { f <- intersect(names(mod_ref), names(m)); ari(mod_ref[f], m[f]) }
  coassign <- function(m) sapply(sort(unique(mod_ref)), function(k) {
    f <- intersect(names(mod_ref)[mod_ref == k], names(m)); tb <- table(m[f])
    sum(tb * (tb - 1) / 2) / (length(f) * (length(f) - 1) / 2)
  })
  q_ratio <- function(g, seeds = 1:10) {   # modularity of reference partition / best found
    f <- V(g)$name; memb <- ifelse(f %in% names(mod_ref), mod_ref[f], NA)
    memb[is.na(memb)] <- max(mod_ref) + seq_len(sum(is.na(memb)))   # unassigned nodes as singletons
    q_best <- max(sapply(seeds, function(s) { set.seed(s); modularity(cluster_louvain(g, weights = E(g)$weight)) }))
    modularity(g, memb, weights = E(g)$weight) / q_best
  }
  coassign_names <- paste0("coassign_M", sort(unique(mod_ref)))

  # baseline: Louvain run-to-run variability on the full network
  g_full <- build_graph(edges_sig)
  baseline <- bind_rows(lapply(seq_len(N_RUNS_FULL), function(s) {
    m <- louvain_run(g_full, s)
    tibble(run = s, ARI = ari_ref(m), !!!setNames(as.list(coassign(m)), coassign_names))
  }))
  baseline$Q_ratio_full <- q_ratio(g_full)

  # (c) LOIO
  loio_edges <- list(); loio_part <- list()
  for (inf in levels(infant)) {
    keep <- infant != inf
    cd_l <- cor_from_ranks(ranks_of(prot_mat[keep, ]), ranks_of(lip_mat[keep, ]),
                           droplevels(infant[keep]), tp[keep], ga[keep], "pooled")
    loio_edges[[inf]] <- compare_edges(cd_l, paste0("LOIO_", inf))
    g_l  <- build_graph(cd_l %>% filter(FDR < FDR_EDGE, abs(rho) >= RHO_FOR_NETWORK))
    runs <- t(sapply(seq_len(N_RUNS_LOIO), function(s) { m <- louvain_run(g_l, s); c(ARI = ari_ref(m), coassign(m)) }))
    loio_part[[inf]] <- tibble(excluded_infant = inf, n_samples = sum(keep),
                               ARI_median = median(runs[, 1]),
                               !!!setNames(as.list(apply(runs[, -1, drop = FALSE], 2, median)), coassign_names),
                               Q_ratio = q_ratio(g_l))
  }
  loio_edges <- bind_rows(loio_edges)
  loio_part  <- bind_rows(loio_part)

  edge_tbl <- bind_rows(
    edge_tbl,
    loio_edges %>% summarise(analysis = "LOIO_min_over_infants", n_edges = min(n_edges),
                             pct_same_sign = min(pct_same_sign), pct_FDR_same_sign = min(pct_FDR_same_sign),
                             pct_nominal_same_sign = min(pct_nominal_same_sign),
                             median_abs_rho_ref = first(median_abs_rho_ref), median_abs_rho = min(median_abs_rho)),
    loio_edges
  )

  write_tsv(edge_tbl,  file.path(rob_dir, "network_robustness_edges_summary.tsv"))
  write_tsv(loio_part, file.path(rob_dir, "network_robustness_loio_partition.tsv"))
  write_tsv(baseline,  file.path(rob_dir, "network_robustness_louvain_baseline.tsv"))

  print(edge_tbl %>% filter(!str_detect(analysis, "^LOIO_P")), width = Inf)
  cat("\nLouvain run-to-run (full data), ARI vs reference: ",
      round(min(baseline$ARI), 2), "-", round(max(baseline$ARI), 2),
      " | Q reference / best: ", round(baseline$Q_ratio_full[1], 3), "\n", sep = "")
  cat("LOIO, ARI vs reference: ", round(min(loio_part$ARI_median), 2), "-", round(max(loio_part$ARI_median), 2),
      " | Q reference / best: ", round(min(loio_part$Q_ratio), 2), "-", round(max(loio_part$Q_ratio), 2), "\n", sep = "")
  cat("LOIO, median co-assignment per module:\n")
  print(round(apply(loio_part %>% select(starts_with("coassign")), 2, median), 2))
  cat("Robustness outputs saved in: ", rob_dir, "\n", sep = "")
}
