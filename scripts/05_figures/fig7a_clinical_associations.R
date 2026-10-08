# ==============================================================================
# Project: Preterm EV Multi-Omics – Integration (proteomics + lipidomics)
# Figure:  Fig. 7 – Clinical associations of cross-omics modules
# Script:  fig7a_clinical_associations.R
# Author:  Nicole Grinovero
# Date:    Jan 2026
#
# Aim:
#   Generate Fig. 7 clinical integration outputs:
#     (A) Heatmap of CAMERA module–outcome associations (proteins vs lipids),
#         visualized as signed −log10(P) with significance marks from FDR_global.
#     (B) Cytoscape-ready network tables (nodes/edges) only.
#         Styling is done manually in Cytoscape.
#
# Inputs:
#   - outputs/tables/preterm/integration/clinical_camera/modules_camera_all_outcomes_multiomics.tsv
#   - outputs/tables/preterm/integration/clinical_camera/drivers_top2_<Layer>_<Outcome>.tsv
#
# Outputs:
#   - outputs/figures/fig7_clinical_associations/
#       fig7_clinical_associations_camera_heatmap.svg
#   - outputs/tables/preterm/integration/clinical_camera/network_fig7_clinical_associations/
#       fig7_clinical_associations_{nodes,edges}_cytoscape.tsv
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(svglite)
  library(grid)      
})

# ------------------------------------------------------------------------------
# 0) Configuration
# ------------------------------------------------------------------------------

# ---- Figure scope ----
FIG_SCOPE <- "all"  # "brain" or "all"

OUTCOMES_ALL <- c(
  "any_brain_lesion",
  "hemorrhagic_lesion",
  "sepsis",
  "mechanical_ventilation"
)

OUTCOMES_BRAIN <- c(
  "any_brain_lesion",
  "hemorrhagic_lesion"
)

OUTCOMES <- if (FIG_SCOPE == "brain") OUTCOMES_BRAIN else OUTCOMES_ALL
LAYERS   <- c("Protein", "Lipid")

# ---- Pretty outcome labels (plots only) ----
outcome_labels <- c(
  "any_brain_lesion"       = "Any brain lesion",
  "hemorrhagic_lesion"     = "Hemorrhagic lesion",
  "sepsis"                 = "Sepsis",
  "mechanical_ventilation" = "Mechanical ventilation"
)

# ---- Significance (heatmap marks) ----
FDR_STAR <- 0.05     # * if FDR_global < 0.05
FDR_DOT  <- 0.10     # · if FDR_global < 0.10
USE_DOT  <- FALSE    # show · for suggestive signals

# ---- Network selection (Outcome→Module edges) ----
FDR_CUTOFF_NET <- 0.05
FALLBACK_TOPK  <- 0      # 0 = do not force edges if none pass

# ---- Drivers shown per selected module (readability) ----
TOP_DRIVERS_PROT <- 5
TOP_DRIVERS_LIP  <- 5

# ---- Palette ----
COL_POS <- "#B33A3A"
COL_NEG <- "#3B6FB6"

# driver node base colors (optional; for Cytoscape initial styling)
COL_PROT <- "#b9d5e1"
COL_LIP  <- "#f3cfad"

# module colours (stored in the node table as base_color)
MOD_LEVELS <- c("M1","M2","M3","M4","M5")
MOD_COLS <- c(
  "M1" = "#C8D7C1",
  "M2" = "#F1D7C2",
  "M3" = "#C9BEDD",
  "M4" = "#D7D7D7",
  "M5" = "#BFD3DE"
)

# ---- Saving helpers ----
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

# ---- Paths ----
in_master <- here(
  "outputs","tables","preterm","integration","clinical_camera",
  "modules_camera_all_outcomes_multiomics.tsv"
)
drivers_dir <- here("outputs","tables","preterm","integration","clinical_camera")

fig_dir <- here("outputs","figures","fig7_clinical_associations")
net_dir <- here(
  "outputs","tables","preterm","integration","clinical_camera",
  "network_fig7_clinical_associations"
)

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(net_dir, recursive = TRUE, showWarnings = FALSE)

out_heat_svg <- file.path(fig_dir, "fig7_clinical_associations_camera_heatmap.svg")

out_nodes <- file.path(net_dir, "fig7_clinical_associations_nodes_cytoscape.tsv")
out_edges <- file.path(net_dir, "fig7_clinical_associations_edges_cytoscape.tsv")

stopifnot(file.exists(in_master))

# ------------------------------------------------------------------------------
# 1) Load MASTER camera table
# ------------------------------------------------------------------------------

cam <- read_tsv(in_master, show_col_types = FALSE) %>%
  mutate(
    Outcome     = as.character(Outcome),
    Layer       = as.character(Layer),
    ModuleLabel = as.character(ModuleLabel),
    N           = suppressWarnings(as.integer(N)),
    Direction   = as.character(Direction),
    PValue      = suppressWarnings(as.numeric(PValue)),
    FDR_global  = suppressWarnings(as.numeric(FDR_global)),
    signed_score = ifelse(Direction == "Up", 1, -1) * (-log10(PValue))
  ) %>%
  filter(Outcome %in% OUTCOMES, Layer %in% LAYERS)

# module ordering
mod_levels <- cam %>%
  distinct(ModuleLabel) %>%
  pull(ModuleLabel) %>%
  as.character()

if (all(MOD_LEVELS %in% mod_levels)) {
  mod_levels <- MOD_LEVELS
} else {
  mod_levels <- sort(mod_levels)
}

cam <- cam %>%
  mutate(
    Outcome     = factor(Outcome, levels = OUTCOMES),
    Layer       = factor(Layer, levels = c("Protein","Lipid")),
    ModuleLabel = factor(ModuleLabel, levels = mod_levels)
  )

# ------------------------------------------------------------------------------
# 2) Panel A — Heatmap (signed −log10(P) + FDR marks)
# ------------------------------------------------------------------------------

heat_df <- cam %>%
  mutate(
    mark = case_when(
      !is.na(FDR_global) & FDR_global < FDR_STAR ~ "*",
      USE_DOT & !is.na(FDR_global) & FDR_global < FDR_DOT ~ "·",
      TRUE ~ ""
    ),
    signed_score = ifelse(is.finite(signed_score), signed_score, NA_real_)
  )

lim <- max(abs(heat_df$signed_score), na.rm = TRUE)
lim <- ifelse(is.finite(lim) && lim > 0, lim, 1)

p_heat <- ggplot(heat_df, aes(x = Outcome, y = ModuleLabel, fill = signed_score)) +
  geom_tile(color = "black", linewidth = 0.35) +               # black border
  geom_text(aes(label = mark), size = 6, color = "grey10") +
  facet_wrap(~ Layer, ncol = 1) +
  scale_x_discrete(labels = outcome_labels) +
  scale_fill_gradient2(
    low = COL_NEG, mid = "white", high = COL_POS,
    midpoint = 0,
    limits   = c(-lim, lim),
    name     = expression("Signed -log"[10]*"(p-value)")
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid   = element_blank(),
    strip.text   = element_text(size = 12, face = "bold", color = "grey10"),
    axis.text.x  = element_text(angle = 35, hjust = 1, vjust = 1, size = 10, color = "grey10"),
    axis.text.y  = element_text(size = 10, color = "grey10"),
    legend.position = "right",
    legend.title = element_text(color = "grey10"),
    legend.text  = element_text(color = "grey10")
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "left",                             # title on left
      title.theme    = element_text(angle = 90, vjust = 0.5),# vertical title
      barwidth       = unit(0.5, "cm"),
      barheight      = unit(3.0, "cm")
    )
  )

print(p_heat)
save_svg(p_heat, out_heat_svg, width = 3.2, height = 4.8)

# ------------------------------------------------------------------------------
# 3) Panel B — Cytoscape network tables (nodes/edges only)
# ------------------------------------------------------------------------------

# --- select Outcome→Module edges for the network ---
sel_edges <- cam %>%
  group_by(Outcome, Layer) %>%
  group_modify(~{
    df <- .x
    sig <- df %>% filter(!is.na(FDR_global), FDR_global < FDR_CUTOFF_NET)
    if (nrow(sig) > 0) sig else df %>% arrange(PValue) %>% slice_head(n = FALLBACK_TOPK)
  }) %>%
  ungroup() %>%
  arrange(Outcome, Layer, PValue) %>%
  mutate(
    OutcomeNode = paste0("OUT:", as.character(Outcome)),
    ModuleNode  = paste0(as.character(Layer), ":", as.character(ModuleLabel)),
    edge_group  = case_when(
      Direction == "Up"   ~ "Outcome→Module:Up",
      Direction == "Down" ~ "Outcome→Module:Down",
      TRUE ~ "Outcome→Module:NA"
    ),
    edge_weight = -log10(PValue)
  )

# --- read drivers (must contain beta, p, FDR; GeneName optional for proteins) ---
read_drivers <- function(layer, outcome) {
  f <- file.path(drivers_dir, paste0("drivers_top2_", layer, "_", outcome, ".tsv"))
  if (!file.exists(f)) return(tibble())
  read_tsv(f, show_col_types = FALSE) %>%
    mutate(
      Layer       = as.character(Layer),
      Outcome     = as.character(Outcome),
      ModuleLabel = as.character(ModuleLabel),
      beta        = suppressWarnings(as.numeric(beta)),
      p           = suppressWarnings(as.numeric(p)),
      FDR         = suppressWarnings(as.numeric(FDR))
    )
}

drivers_all <- map_dfr(as.character(OUTCOMES), \(oc) {
  bind_rows(
    read_drivers("Protein", oc),
    read_drivers("Lipid",   oc)
  )
})

# keep only drivers for selected modules in sel_edges
sel_mods <- sel_edges %>%
  distinct(
    Outcome = as.character(Outcome),
    Layer = as.character(Layer),
    ModuleLabel = as.character(ModuleLabel)
  ) %>%
  mutate(ModuleNode = paste0(Layer, ":", ModuleLabel))

drivers_all <- drivers_all %>%
  mutate(
    Outcome     = as.character(Outcome),
    Layer       = as.character(Layer),
    ModuleLabel = as.character(ModuleLabel),
    ModuleNode  = paste0(Layer, ":", ModuleLabel)
  ) %>%
  inner_join(sel_mods, by = c("Outcome","Layer","ModuleLabel","ModuleNode"))

# label + stable node IDs (proteins: prefer GeneName if present)
drivers_all <- drivers_all %>%
  mutate(
    DriverLabel = case_when(
      Layer == "Protein" &
        "GeneName" %in% names(.) &
        !is.na(GeneName) & GeneName != "" ~ as.character(GeneName),
      TRUE ~ as.character(Feature)
    ),
    FeatureNode = ifelse(Layer == "Protein",
                         paste0("P:", DriverLabel),
                         paste0("L:", DriverLabel))
  )

# limit drivers per (Outcome×Layer×Module) for readability
drivers_sel <- bind_rows(
  drivers_all %>%
    filter(Layer == "Protein") %>%
    group_by(Outcome, Layer, ModuleNode) %>%
    arrange(p, desc(abs(beta)), .by_group = TRUE) %>%
    slice_head(n = TOP_DRIVERS_PROT) %>%
    ungroup(),
  drivers_all %>%
    filter(Layer == "Lipid") %>%
    group_by(Outcome, Layer, ModuleNode) %>%
    arrange(p, desc(abs(beta)), .by_group = TRUE) %>%
    slice_head(n = TOP_DRIVERS_LIP) %>%
    ungroup()
)

# --- edges table ---
edges_outcome_module <- sel_edges %>%
  transmute(
    source      = OutcomeNode,
    target      = ModuleNode,
    interaction = "assoc",
    edge_type   = "Outcome→Module",
    edge_group  = edge_group,
    Outcome     = as.character(Outcome),
    Layer       = as.character(Layer),
    ModuleLabel = as.character(ModuleLabel),
    Direction   = as.character(Direction),
    PValue      = as.numeric(PValue),
    FDR_global  = as.numeric(FDR_global),
    signed_score = as.numeric(signed_score),
    weight      = as.numeric(edge_weight)
  )

edges_module_driver <- drivers_sel %>%
  transmute(
    source      = ModuleNode,
    target      = FeatureNode,
    interaction = "contains",
    edge_type   = "Module→Driver",
    edge_group  = "Module→Driver",
    Outcome     = as.character(Outcome),
    Layer       = as.character(Layer),
    ModuleLabel = as.character(ModuleLabel),
    Feature     = as.character(Feature),
    DriverLabel = as.character(DriverLabel),
    beta        = as.numeric(beta),
    p           = as.numeric(p),
    FDR         = as.numeric(FDR),
    weight      = abs(beta)
  )

edges_cyto <- bind_rows(edges_outcome_module, edges_module_driver)

# --- nodes table ---
nodes_outcomes <- edges_outcome_module %>%
  distinct(id = source, Outcome) %>%
  transmute(
    id,
    label = Outcome,
    node_type = "Outcome",
    layer = NA_character_,
    module = NA_character_,
    feature_layer = NA_character_,
    beta = NA_real_,
    base_color = "#444444",
    base_shape = "DIAMOND"
  )

nodes_modules <- edges_outcome_module %>%
  distinct(id = target, Layer, ModuleLabel) %>%
  transmute(
    id,
    label = ModuleLabel,
    node_type = "Module",
    layer = Layer,
    module = ModuleLabel,
    feature_layer = NA_character_,
    beta = NA_real_,
    base_color = ifelse(ModuleLabel %in% names(MOD_COLS), MOD_COLS[ModuleLabel], "#BBBBBB"),
    base_shape = "ELLIPSE"
  )

nodes_drivers <- edges_module_driver %>%
  distinct(id = target, Layer, DriverLabel, beta) %>%
  transmute(
    id,
    label = DriverLabel,
    node_type = ifelse(Layer == "Protein", "Protein", "Lipid"),
    layer = NA_character_,
    module = NA_character_,
    feature_layer = Layer,
    beta = as.numeric(beta),              # keep numeric in file
    base_color = ifelse(Layer == "Protein", COL_PROT, COL_LIP),
    base_shape = ifelse(Layer == "Protein", "ROUND_RECTANGLE", "RECTANGLE")
  )

nodes_cyto <- bind_rows(nodes_outcomes, nodes_modules, nodes_drivers) %>%
  distinct(id, .keep_all = TRUE)

write_tsv(nodes_cyto, out_nodes)
write_tsv(edges_cyto, out_edges)

cat(
  "\n[OK] Fig. 7 outputs written:\n",
  " - Heatmap: ", out_heat_svg, "\n",
  " - Cytoscape nodes: ", out_nodes, "\n",
  " - Cytoscape edges: ", out_edges, "\n\n",
  "Cytoscape import:\n",
  "  1) File → Import → Network from File...  (", basename(out_edges), ")\n",
  "  2) File → Import → Table → File... (", basename(out_nodes), "; key = id)\n\n",
  "Manual Cytoscape mapping tips:\n",
  "  - Node Fill Color: map 'beta' (continuous) for driver nodes\n",
  "  - Edge Color: map 'edge_group' (discrete)\n",
  "  - Edge Width: map 'weight' (continuous)\n",
  sep = ""
)
