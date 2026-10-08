# ==============================================================================
# Project: Preterm EV Multi-Omics – Proteomics
# Script:  06_prot_enrich_T4vsT0_GOBP_Reactome.R
# Author:  Nicole Grinovero
# Date:    Jan 2026
#
# Aim:
#   Enrichment analysis to assess GO Biological Process and Reactome pathways
#   among proteins changing at term-equivalent age (T4) vs birth (T0),
#   separately for UP- and DOWN-regulated proteins (defined by LMM results).
#
# Input:
#   - outputs/tables/preterm/proteomics/proteins_LMM_time_results.txt
#
# Output:
#   - outputs/tables/preterm/proteomics/enrich_GOBP_Reactome_T4vsT0.txt
#
# Notes:
#   - Universe = all measured proteins in the LMM results (GeneName non-missing)
#   - UP/DOWN defined by FDR + effect size threshold (beta)
#   - Enrichment uses Entrez IDs (SYMBOL -> ENTREZID mapping via org.Hs.eg.db)
#   - QC: script prints mapping and list sizes to help reproducibility/debugging
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
})

# ------------------------------------------------------------------------------
# 0) Configuration
# ------------------------------------------------------------------------------
stats_path <- here("outputs", "tables", "preterm", "proteomics", "proteins_LMM_time_results.txt")

out_dir <- here("outputs", "tables", "preterm", "proteomics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, "enrich_GOBP_Reactome_T4vsT0.txt")

# Define UP/DOWN sets from LMM
BETA_CUTOFF <- log2(1.5)
FDR_CUTOFF  <- 0.05

# Gene set size filters for GO; MIN_GS is also the minimum list size for
# running GO and Reactome enrichment
MIN_GS <- 10
MAX_GS <- 500

# GO term redundancy reduction (semantic similarity)
SIM_CUTOFF <- 0.7

# Top terms retained per database and direction
N_GO <- 3
N_RE <- 3

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------
map_symbols_to_entrez <- function(symbols, label = "set") {
  symbols <- unique(str_squish(symbols))
  symbols <- symbols[!is.na(symbols) & symbols != ""]
  if (length(symbols) == 0) {
    return(list(entrez = character(0), map = tibble(), n_in = 0L, n_mapped = 0L))
  }
  
  conv <- bitr(
    symbols,
    fromType = "SYMBOL", toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  ) %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(ENTREZID, .keep_all = TRUE)
  
  entrez <- unique(conv$ENTREZID)
  
  message("[Mapping] ", label, ": symbols=", length(symbols),
          " | mapped_entrez=", length(entrez),
          " | unmapped_symbols~=", max(0, length(symbols) - nrow(conv)))
  
  list(entrez = entrez, map = conv, n_in = length(symbols), n_mapped = length(entrez))
}

safe_enrichGO <- function(gene, universe) {
  if (length(gene) < MIN_GS) return(NULL)
  
  enrichGO(
    gene          = gene,
    universe      = universe,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff  = FDR_CUTOFF,
    minGSSize     = MIN_GS,
    maxGSSize     = MAX_GS,
    readable      = TRUE
  )
}

safe_simplify_GO <- function(ego) {
  if (is.null(ego)) return(NULL)
  df <- tryCatch(as.data.frame(ego), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(ego)
  simplify(ego, cutoff = SIM_CUTOFF, by = "p.adjust", select_fun = min)
}

safe_enrichReactome <- function(gene, universe) {
  if (length(gene) < MIN_GS) return(NULL)
  
  enrichPathway(
    gene          = gene,
    universe      = universe,
    organism      = "human",
    pAdjustMethod = "BH",
    qvalueCutoff  = FDR_CUTOFF,
    readable      = TRUE
  )
}

pick_top <- function(df, n, source_name, direction) {
  if (is.null(df) || nrow(df) == 0) return(tibble())
  df %>%
    arrange(p.adjust, desc(Count)) %>%
    slice_head(n = min(n, nrow(df))) %>%
    transmute(
      Direction = direction,
      Source    = source_name,
      Term      = str_to_sentence(str_replace_all(Description, "_", " ")),
      FDR       = p.adjust,
      Count     = Count,
      neglog10  = -log10(p.adjust),
      Genes     = geneID
    )
}

# ------------------------------------------------------------------------------
# 1) Load statistics and define universe, UP and DOWN gene sets (SYMBOL)
# ------------------------------------------------------------------------------
stats <- read_delim(stats_path, delim = "\t", show_col_types = FALSE)

stopifnot(all(c("GeneName", "beta_T4_vs_T0", "FDR_T4_T0") %in% names(stats)))

stats2 <- stats %>%
  transmute(
    Gene = str_squish(as.character(GeneName)),
    beta = as.numeric(beta_T4_vs_T0),
    fdr  = as.numeric(FDR_T4_T0)
  ) %>%
  filter(!is.na(Gene), Gene != "", !is.na(beta), !is.na(fdr))

universe_symbols <- unique(stats2$Gene)

up_symbols <- stats2 %>%
  filter(fdr < FDR_CUTOFF, beta >  BETA_CUTOFF) %>%
  pull(Gene) %>%
  unique()

down_symbols <- stats2 %>%
  filter(fdr < FDR_CUTOFF, beta < -BETA_CUTOFF) %>%
  pull(Gene) %>%
  unique()

message("Universe (symbols): ", length(universe_symbols))
message("UP (symbols):       ", length(up_symbols))
message("DOWN (symbols):     ", length(down_symbols))

if (length(up_symbols) == 0 && length(down_symbols) == 0) {
  warning("No UP or DOWN genes passed thresholds. Output will be empty.")
}

# ------------------------------------------------------------------------------
# 2) Convert SYMBOL -> ENTREZID (required by clusterProfiler/ReactomePA)
# ------------------------------------------------------------------------------
conv_univ <- map_symbols_to_entrez(universe_symbols, label = "universe")
conv_up   <- map_symbols_to_entrez(up_symbols,       label = "UP")
conv_dn   <- map_symbols_to_entrez(down_symbols,     label = "DOWN")

universe_entrez <- conv_univ$entrez
up_entrez       <- conv_up$entrez
down_entrez     <- conv_dn$entrez

# ------------------------------------------------------------------------------
# 3) Enrichment: GO Biological Process (with GO redundancy reduction)
# ------------------------------------------------------------------------------
ego_up <- safe_enrichGO(up_entrez, universe_entrez)
ego_dn <- safe_enrichGO(down_entrez, universe_entrez)

ego_up <- safe_simplify_GO(ego_up)
ego_dn <- safe_simplify_GO(ego_dn)

res_go_up <- if (!is.null(ego_up)) as.data.frame(ego_up) %>% filter(p.adjust < FDR_CUTOFF) else tibble()
res_go_dn <- if (!is.null(ego_dn)) as.data.frame(ego_dn) %>% filter(p.adjust < FDR_CUTOFF) else tibble()

# ------------------------------------------------------------------------------
# 4) Enrichment: Reactome pathways
# ------------------------------------------------------------------------------
er_up <- safe_enrichReactome(up_entrez, universe_entrez)
er_dn <- safe_enrichReactome(down_entrez, universe_entrez)

res_re_up <- if (!is.null(er_up)) as.data.frame(er_up) %>% filter(p.adjust < FDR_CUTOFF) else tibble()
res_re_dn <- if (!is.null(er_dn)) as.data.frame(er_dn) %>% filter(p.adjust < FDR_CUTOFF) else tibble()

# ------------------------------------------------------------------------------
# 5) Build a single tidy output table (top terms only)
# ------------------------------------------------------------------------------
enrich_tbl <- bind_rows(
  pick_top(res_go_up, N_GO, "GO:BP",    "Up"),
  pick_top(res_re_up, N_RE, "Reactome", "Up"),
  pick_top(res_go_dn, N_GO, "GO:BP",    "Down"),
  pick_top(res_re_dn, N_RE, "Reactome", "Down")
)

write_tsv(enrich_tbl, out_file)

cat("\n[OK] Enrichment table saved:\n - ", out_file, "\n", sep = "")
cat("[OK] Rows saved: ", nrow(enrich_tbl), "\n", sep = "")
