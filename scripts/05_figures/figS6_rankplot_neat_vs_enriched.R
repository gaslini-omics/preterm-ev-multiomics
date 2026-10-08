# ==============================================================================
# Project: Preterm EV Multi-Omics – Technical validation
# Figure:  Fig. S6 – EV marker enrichment: rank-abundance, neat vs Mag-Net BAMM 10 µL
# Script:  figS6_rankplot_neat_vs_enriched.R
# Author:  Nicole Grinovero
#
# Aim:
#   Rank-abundance plots of the plasma proteome analyzed directly (neat) or after
#   Mag-Net BAMM enrichment from 10 µL of pooled adult plasma, with proteins
#   coloured by MISEV2023 category (Table 3 of the guidelines):
#     - EV markers        = categories 1-2
#     - Non-EV co-isolates = category 3
#   Within each preparation, proteins quantified in >= 3 of 4 technical
#   replicates are ranked by mean intensity (no imputation).
#
# Inputs:
#   - data/techrep/proteomics/20260206_182115_plasma_EV_multiomics_techrep_Report.tsv
#   - data/reference/misev2023_category_assignment.tsv   (GeneName, MISEVcat)
#
# Outputs:
#   - outputs/figures/figS6/figS6_rankplot_neat_vs_enriched.{svg,png}
#   - outputs/figures/figS6/figS6_source_data.tsv
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(ggrepel)
  library(svglite)
})

# ------------------------------------------------------------------------------
# 0) Paths and settings
# ------------------------------------------------------------------------------
report_path <- here("data", "techrep", "proteomics",
                    "20260206_182115_plasma_EV_multiomics_techrep_Report.tsv")
misev_path  <- here("data", "reference", "misev2023_category_assignment.tsv")

out_dir <- here("outputs", "figures", "figS6")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

MIN_REP <- 3   # quantified in >= 3 of 4 replicates

conditions <- c("Neat" = "Neat",
                "MagNet_BAMM_10uL" = "Mag-Net BAMM 10µL")

cat_levels <- c("EV markers (MISEV cat. 1-2)",
                "Non-EV co-isolates (MISEV cat. 3)",
                "Other proteins")
cat_cols   <- c("EV markers (MISEV cat. 1-2)"       = "#A25E6C",
                "Non-EV co-isolates (MISEV cat. 3)" = "#8FA6B5",
                "Other proteins"                    = "#A6A6A6")

# proteins labelled in each panel
label_genes <- list(
  "Neat"              = c("ALB", "APOA1", "APOB", "LAMP1"),
  "Mag-Net BAMM 10µL" = c("ALB", "APOA1", "APOB", "ITGB1", "HLA-G", "ADAM10",
                          "PDCD6IP", "SDC4", "FLOT1", "FLOT2", "CHMP4B", "TSG101")
)

# ------------------------------------------------------------------------------
# 1) Helpers
# ------------------------------------------------------------------------------
# gene symbol of the first Swiss-Prot ("sp") entry of a protein group
extract_id <- function(db_str, id_str, target_db = "sp") {
  if (is.na(db_str) || is.na(id_str)) return(NA_character_)
  db_values <- str_trim(str_split(db_str, ";", simplify = TRUE))
  id_values <- str_trim(str_split(id_str, ";", simplify = TRUE))
  if (length(id_values) == 1) return(id_values[1])
  idx <- which(db_values == target_db)
  if (length(idx) > 0) return(id_values[idx[1]])
  id_values[1]
}

# ------------------------------------------------------------------------------
# 2) Load data
# ------------------------------------------------------------------------------
rep <- read_delim(report_path, delim = "\t", show_col_types = FALSE)
misev <- read_tsv(misev_path, show_col_types = FALSE)

sample_cols <- grep("\\.raw\\.PG\\.Quantity$", names(rep), value = TRUE)

long <- rep %>%
  mutate(GeneName = map2_chr(PG.Database, PG.Genes, extract_id)) %>%
  select(PG.ProteinGroups, GeneName, all_of(sample_cols)) %>%
  pivot_longer(all_of(sample_cols), names_to = "sample", values_to = "Intensity") %>%
  mutate(
    Intensity = suppressWarnings(as.numeric(Intensity)),
    cond_raw  = sample %>%
      str_remove("^\\[\\d+\\]\\s*") %>%
      str_remove("_C\\d+_\\d+\\.raw\\.PG\\.Quantity$")
  ) %>%
  filter(cond_raw %in% names(conditions), !is.na(GeneName), GeneName != "") %>%
  mutate(Condition = unname(conditions[cond_raw]))

# ------------------------------------------------------------------------------
# 3) Per-condition abundance, detection filter, one protein group per gene, rank
# ------------------------------------------------------------------------------
rank_tbl <- long %>%
  group_by(Condition, GeneName, PG.ProteinGroups) %>%
  summarise(
    n_det = sum(is.finite(Intensity) & Intensity > 0),
    mu    = mean(Intensity[is.finite(Intensity) & Intensity > 0]),
    .groups = "drop"
  ) %>%
  filter(n_det >= MIN_REP) %>%
  arrange(Condition, GeneName, desc(n_det), desc(mu)) %>%
  group_by(Condition, GeneName) %>%
  slice(1) %>%
  ungroup() %>%
  group_by(Condition) %>%
  arrange(desc(mu), .by_group = TRUE) %>%
  mutate(Rank = row_number()) %>%
  ungroup() %>%
  mutate(log2I = log2(mu + 1)) %>%   # pseudocount of 1
  left_join(misev, by = "GeneName") %>%
  mutate(
    Category = case_when(
      MISEVcat %in% c(1, 2) ~ "EV markers (MISEV cat. 1-2)",
      MISEVcat == 3         ~ "Non-EV co-isolates (MISEV cat. 3)",
      TRUE                  ~ "Other proteins"
    ),
    Category  = factor(Category, levels = cat_levels),
    Condition = factor(Condition, levels = unname(conditions))
  ) %>%
  select(Condition, GeneName, n_det, mu, log2I, Rank, Category, MISEVcat) %>%
  arrange(Condition, Rank)

write_tsv(rank_tbl, file.path(out_dir, "figS6_source_data.tsv"))

cat("ALB rank:\n")
print(rank_tbl %>% filter(GeneName == "ALB") %>% select(Condition, Rank))

# ------------------------------------------------------------------------------
# 4) Plot
# ------------------------------------------------------------------------------
counts <- rank_tbl %>%
  group_by(Condition) %>%
  summarise(n_prot = n(),
            n_ev   = sum(Category == cat_levels[1]), .groups = "drop") %>%
  mutate(lab_prot = paste(n_prot, "proteins"),
         lab_ev   = paste(n_ev, "EV markers"))

labels_df <- rank_tbl %>%
  filter(map2_lgl(as.character(Condition), GeneName,
                  ~ .y %in% label_genes[[.x]]))

bg <- rank_tbl %>% filter(Category == "Other proteins")
fg <- rank_tbl %>% filter(Category != "Other proteins")

p <- ggplot(rank_tbl, aes(Rank, log2I)) +
  geom_point(data = bg, colour = cat_cols[["Other proteins"]],
             size = 1.1, alpha = 0.18) +
  geom_point(data = fg, aes(fill = Category), shape = 21, colour = "grey20",
             stroke = 0.3, size = 2, alpha = 0.95) +
  geom_text_repel(data = labels_df, aes(label = GeneName), size = 3.2,
                  segment.colour = "grey40", segment.size = 0.3,
                  min.segment.length = 0, box.padding = 0.4,
                  max.overlaps = Inf, seed = 1) +
  geom_text(data = counts, aes(x = Inf, y = Inf, label = lab_prot),
            hjust = 1.1, vjust = 1.6, size = 3.3, colour = "#595959",
            fontface = "bold", inherit.aes = FALSE) +
  geom_text(data = counts, aes(x = Inf, y = Inf, label = lab_ev),
            hjust = 1.1, vjust = 3.4, size = 3.3, colour = "#9B5968",
            inherit.aes = FALSE) +
  facet_wrap(~ Condition, nrow = 1) +
  scale_fill_manual(values = cat_cols[cat_levels[1:2]], name = NULL,
                    labels = c(expression("EV markers " * italic("(MISEV cat. 1-2)")),
                               expression("Non-EV co-isolates " * italic("(MISEV cat. 3)")))) +
  scale_x_continuous(limits = c(0, 1000)) +
  labs(x = "Protein rank", y = expression(log[2] * "(mean intensity)")) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid       = element_blank(),
    strip.background = element_rect(fill = "#EBEBEB", colour = "grey20"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "top",
    text             = element_text(family = "Arial")
  )

ggsave(file.path(out_dir, "figS6_rankplot_neat_vs_enriched.svg"), p,
       width = 7, height = 4.2, device = svglite::svglite)
ggsave(file.path(out_dir, "figS6_rankplot_neat_vs_enriched.png"), p,
       width = 7, height = 4.2, dpi = 300, bg = "white")

message("Fig. S6 written to: ", out_dir)
