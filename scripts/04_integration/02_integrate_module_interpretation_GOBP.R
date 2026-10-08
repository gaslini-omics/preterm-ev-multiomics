# ==============================================================================
# Project: Preterm EV Multi-Omics – Integration (proteomics + lipidomics)
# Script:  02_integrate_module_interpretation_GOBP.R
# Author:  Nicole Grinovero
# Date:    Jan 2026
#
# Aim:
#   GO:BP enrichment of cross-omics modules (protein members only), with an
#   option to run enrichment either:
#     (A) using a measured-proteome universe (background), or
#     (B) without universe (default GO background),
#   consistently for ALL modules.
#
#   Outputs are tidy and include per-term associated genes (symbols) for plotting.
#
# Inputs:
#   - outputs/tables/preterm/integration/integration_modules_membership_annotated.tsv
#   - outputs/tables/preterm/proteomics/proteins_LMM_time_results.txt
#
# Outputs:
#   - outputs/tables/preterm/integration/integration_modules_enrichment_GOBP.tsv
#   - outputs/tables/preterm/integration/integration_modules_interpretation_summary.tsv
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# ------------------------------------------------------------------------------
# 0) Configuration
# ------------------------------------------------------------------------------
USE_BACKGROUND <- FALSE    # TRUE = measured-proteome universe; FALSE = default GO background
FILTER_UNIVERSE_TO_BP <- FALSE  # recommended if USE_BACKGROUND=TRUE

FDR_CUTOFF <- 0.05
MIN_GS     <- 10
MAX_GS     <- 1000
SIM_CUTOFF <- 0.7
TOP_GO     <- 5
MIN_PROT_MAPPED <- 8

PROT_GENE <- "GeneName"
LIP_CLASS <- "Class"
LIP_CAT   <- "Category"

modules_path   <- here("outputs","tables","preterm","integration","integration_modules_membership_annotated.tsv")
prot_stats_path <- here("outputs","tables","preterm","proteomics","proteins_LMM_time_results.txt")

out_dir  <- here("outputs","tables","preterm","integration")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_main <- file.path(out_dir, "integration_modules_enrichment_GOBP.tsv")
out_sum  <- file.path(out_dir, "integration_modules_interpretation_summary.tsv")

# ------------------------------------------------------------------------------
# 1) Helpers (minimal)
# ------------------------------------------------------------------------------
map_symbols_to_entrez <- function(symbols) {
  symbols <- unique(str_squish(as.character(symbols)))
  symbols <- symbols[!is.na(symbols) & symbols != ""]
  if (!length(symbols)) return(character(0))
  
  conv <- suppressMessages(
    bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  )
  unique(conv$ENTREZID)
}

simplify_go <- function(ego) {
  if (is.null(ego)) return(NULL)
  simplify(ego, cutoff = SIM_CUTOFF, by = "p.adjust", select_fun = min)
}

collapse_top <- function(x, n = 3) {
  x <- str_squish(as.character(x))
  x <- x[!is.na(x) & x != ""]
  if (!length(x)) return(NA_character_)
  paste0(head(names(sort(table(x), decreasing = TRUE)), n), collapse = "; ")
}

# ------------------------------------------------------------------------------
# 2) Load inputs
# ------------------------------------------------------------------------------
mods <- read_tsv(modules_path, show_col_types = FALSE) %>%
  mutate(
    ModuleLabel = as.character(ModuleLabel),
    Module      = as.integer(Module),
    Type        = as.character(Type),
    Gene        = str_squish(as.character(.data[[PROT_GENE]]))
  )

prot_stats <- read_delim(prot_stats_path, delim = "\t", show_col_types = FALSE)

# ------------------------------------------------------------------------------
# 3) Universe (optional background)
# ------------------------------------------------------------------------------
universe_entrez <- map_symbols_to_entrez(prot_stats[[PROT_GENE]])

if (USE_BACKGROUND && FILTER_UNIVERSE_TO_BP) {
  bp_map <- suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db,
    keys = universe_entrez,
    keytype = "ENTREZID",
    columns = c("GO","ONTOLOGY")
  ))
  universe_entrez <- bp_map %>%
    filter(ONTOLOGY == "BP", !is.na(GO)) %>%
    pull(ENTREZID) %>%
    unique()
}

# ------------------------------------------------------------------------------
# 4) Module composition + lipid summary + hubs
# ------------------------------------------------------------------------------
mod_comp <- mods %>%
  count(ModuleLabel, Module, Type, name = "n") %>%
  pivot_wider(names_from = Type, values_from = n, values_fill = 0) %>%
  transmute(
    ModuleLabel, Module,
    n_Protein = ifelse(is.na(Protein), 0, Protein),
    n_Lipid   = ifelse(is.na(Lipid),   0, Lipid),
    n_Total   = n_Protein + n_Lipid
  ) %>%
  arrange(Module)

lip_top <- mods %>%
  filter(Type == "Lipid") %>%
  mutate(
    LipClass = if (LIP_CLASS %in% names(mods)) as.character(.data[[LIP_CLASS]]) else NA_character_,
    LipCat   = if (LIP_CAT   %in% names(mods)) as.character(.data[[LIP_CAT]])   else NA_character_
  ) %>%
  group_by(ModuleLabel, Module) %>%
  summarise(
    TopClass    = collapse_top(LipClass, n = 3),
    TopCategory = collapse_top(LipCat,   n = 3),
    .groups = "drop"
  )

hubs_top <- mods %>%
  mutate(Degree = suppressWarnings(as.numeric(.data$Degree))) %>%
  group_by(ModuleLabel, Module) %>%
  arrange(desc(Degree)) %>%
  summarise(TopHubs = paste0(head(Feature, 5), collapse = "; "), .groups = "drop")

# ------------------------------------------------------------------------------
# 5) GO:BP enrichment per module (proteins only)
# ------------------------------------------------------------------------------
prot_members <- mods %>%
  filter(Type == "Protein", !is.na(Gene), Gene != "") %>%
  dplyr::select(ModuleLabel, Module, Gene) %>%
  distinct()

modules_to_test <- mod_comp$Module

enrich_tbl <- map_dfr(modules_to_test, function(m) {
  
  lab <- paste0("M", m)
  
  genes_m  <- prot_members %>% filter(Module == m) %>% pull(Gene) %>% unique()
  entrez_m <- map_symbols_to_entrez(genes_m)
  
  if (length(entrez_m) < MIN_PROT_MAPPED) {
    return(tibble(
      ModuleLabel = lab, Module = m,
      GO_ID = NA_character_, GO_Term = NA_character_,
      GO_FDR = NA_real_, GO_p = NA_real_, GO_Count = NA_integer_,
      GO_neglog10 = NA_real_, is_FDR_sig = NA,
      TermGenes_Entrez = NA_character_, TermGenes_Symbol = NA_character_,
      n_symbols = length(genes_m), n_entrez = length(entrez_m)
    ))
  }
  
  ego <- if (USE_BACKGROUND) {
    enrichGO(
      gene = entrez_m, universe = universe_entrez,
      OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
      pAdjustMethod = "BH", qvalueCutoff = 1,
      minGSSize = MIN_GS, maxGSSize = MAX_GS,
      readable = FALSE
    )
  } else {
    enrichGO(
      gene = entrez_m,
      OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
      pAdjustMethod = "BH", qvalueCutoff = 1,
      minGSSize = MIN_GS, maxGSSize = MAX_GS,
      readable = FALSE
    )
  }
  
  ego <- simplify_go(ego)
  res <- as.data.frame(ego)
  
  res_top <- res %>%
    arrange(p.adjust, desc(Count)) %>%
    slice_head(n = min(TOP_GO, nrow(res)))
  
  gene_ids_entrez <- unique(unlist(strsplit(res_top$geneID, "/", fixed = TRUE)))
  gene_ids_entrez <- gene_ids_entrez[gene_ids_entrez != ""]
  
  sym_map <- suppressMessages(
    bitr(gene_ids_entrez, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
  ) %>%
    filter(!is.na(SYMBOL)) %>%
    distinct(ENTREZID, .keep_all = TRUE)
  
  map_term_genes_to_symbols <- function(geneID_str) {
    ids <- unlist(strsplit(geneID_str, "/", fixed = TRUE))
    ids <- ids[ids != ""]
    syms <- sym_map$SYMBOL[match(ids, sym_map$ENTREZID)]
    syms <- unique(syms[!is.na(syms) & syms != ""])
    paste0(syms, collapse = ";")
  }
  
  res_top %>%
    transmute(
      ModuleLabel = lab,
      Module      = m,
      GO_ID       = ID,
      GO_Term     = str_to_sentence(str_replace_all(Description, "_", " ")),
      GO_FDR      = p.adjust,
      GO_p        = pvalue,
      GO_Count    = Count,
      GO_neglog10 = -log10(p.adjust),
      is_FDR_sig  = (p.adjust < FDR_CUTOFF),
      TermGenes_Entrez  = geneID,
      TermGenes_Symbol  = vapply(geneID, map_term_genes_to_symbols, character(1)),
      n_symbols   = length(genes_m),
      n_entrez    = length(entrez_m)
    )
})

# ------------------------------------------------------------------------------
# 6) Save MAIN + summary
# ------------------------------------------------------------------------------
main_tbl <- enrich_tbl %>%
  left_join(mod_comp, by = c("ModuleLabel","Module")) %>%
  left_join(lip_top,  by = c("ModuleLabel","Module")) %>%
  left_join(hubs_top, by = c("ModuleLabel","Module")) %>%
  arrange(Module, GO_FDR)

write_tsv(main_tbl, out_main)

sum_tbl <- main_tbl %>%
  group_by(ModuleLabel, Module, n_Total, n_Protein, n_Lipid, TopClass, TopCategory, TopHubs) %>%
  summarise(
    TopGO = {
      x <- GO_Term[!is.na(GO_Term)]
      if (!length(x)) NA_character_ else paste0(head(x, 3), collapse = "; ")
    },
    Any_FDR_sig = any(is_FDR_sig %in% TRUE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Module)

write_tsv(sum_tbl, out_sum)

cat("\n[OK] Saved:\n - ", out_main, "\n - ", out_sum, "\n", sep = "")
cat("[OK] USE_BACKGROUND = ", USE_BACKGROUND, "\n", sep = "")
