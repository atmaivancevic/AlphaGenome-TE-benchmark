#!/usr/bin/env Rscript
#
# Per-track AlphaGenome raw_score + quantile_score for the four headline Fig 3
# candidates, exposing within-protocol RNA-seq variability (stranded total /
# stranded polyA / unstranded polyA) — most agree, but small inserts (e.g. the
# NR1H3 53 bp SVA) can show stranded-vs-unstranded disagreement. Outputs a
# multi-panel vector PDF + per-variant PDFs.
#
# Example usage:
# Rscript scripts/fig3_polymorphic_TE_eQTLs/plot_supp_fig_AG_per_track.R

suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(readr); library(ggplot2); library(patchwork)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--ag_dir",  default = "results/AG_LFC_polymorphic_TE"),
  make_option("--out_dir", default = "figures/supp_fig_AG_per_track")
)))
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# Four headline candidates — locked from tasks/stage4_task4.3_fig3_eqtl_pipeline_locked.md
candidates <- tibble::tribble(
  ~variant_id,         ~family,    ~gene,      ~order,
  "SvimAsm00107100",   "Alu",      "HSD17B12",  1L,
  "SvimAsm00022857",   "L1",       "NEMP2",     2L,
  "SvimAsm00107233",   "SVA",      "NR1H3",     3L,
  "SvimAsm00060017",   "LTR5_Hs",  "HLA-DQA2",  4L
)

# Two-track decomposition (Liu 2026): stranded total + stranded polyA RNA-seq;
# the unstranded polyA track (strand == '.') is dropped.
TRACK_COLORS <- c(
  "total RNA-seq (stranded)"     = "#3A7CA5",   # blue
  "polyA RNA-seq (stranded)"     = "#5BA75A"    # green
)

label_track <- function(track_name, track_strand) {
  is_total <- grepl("total RNA-seq", track_name)
  ifelse(is_total,
         "total RNA-seq (stranded)",
         "polyA RNA-seq (stranded)")
}

build_one_variant_panel <- function(vid, fam, gene) {
  csv <- file.path(opt$ag_dir, sprintf("%s_GM12878.csv", vid))
  ag <- read_csv(csv, show_col_types = FALSE)
  rows <- ag %>%
    filter(gene_name == gene, track_strand %in% c("+", "-")) %>%
    mutate(track_label = label_track(track_name, track_strand)) %>%
    select(track_label, raw_score, quantile_score)
  if (nrow(rows) == 0) {
    stop(sprintf("No rows for %s/%s in %s", vid, gene, csv))
  }
  rows$track_label <- factor(rows$track_label, levels = names(TRACK_COLORS))

  # Long format for faceting raw vs quantile
  long <- bind_rows(
    rows %>% transmute(track_label, metric = "raw score", value = raw_score),
    rows %>% transmute(track_label, metric = "quantile score", value = quantile_score)
  )
  long$metric <- factor(long$metric, levels = c("raw score", "quantile score"))

  title <- sprintf("%s - %s (%s)", vid, gene, fam)

  p <- ggplot(long, aes(x = track_label, y = value, fill = track_label)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, color = "grey30", linewidth = 0.3) +
    facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = TRACK_COLORS, drop = FALSE, guide = "none") +
    labs(title = title, x = NULL, y = NULL) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11),
          axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
          strip.text = element_text(face = "bold", size = 9.5),
          panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3))

  # Save per-variant PDF (cairo for vector text)
  ggsave(file.path(opt$out_dir, sprintf("%s_%s.pdf", vid,
                                        gsub("[^A-Za-z0-9]", "_", gene))),
         p, width = 6.5, height = 3.0, device = "pdf", useDingbats = FALSE)

  p
}

panels <- candidates %>% arrange(order) %>%
  rowwise() %>% mutate(panel = list(build_one_variant_panel(variant_id, family, gene))) %>%
  pull(panel)

# Composite: 4 rows × 1 col (each row holds the 2-facet (raw, quantile) panel for one variant)
composite <- wrap_plots(panels, ncol = 1) +
  plot_annotation(
    title = "Supplementary Figure SX. Per-track AlphaGenome RNA-seq predictions for the four headlined polymorphic-TE variants in GM12878.",
    subtitle = "Stranded total (rRNA-depleted) and stranded polyA RNA-seq tracks (Liu et al 2026 convention). Per-track decomposition of the mean values reported in Supplementary Tables 4 and 5.",
    theme = theme(plot.title = element_text(size = 11, face = "bold"),
                  plot.subtitle = element_text(size = 9, face = "italic", color = "grey25"))
  )

out_pdf <- file.path(opt$out_dir, "supp_fig_AG_per_track.pdf")
ggsave(out_pdf, composite, width = 7.5, height = 11.5, device = "pdf", useDingbats = FALSE)
cat(sprintf("Wrote %s\n", out_pdf))

# Also dump the underlying tabular data alongside the PDF for transparency.
tabular <- candidates %>% arrange(order) %>%
  rowwise() %>% do({
    vid <- .$variant_id; fam <- .$family; gene <- .$gene
    ag <- read_csv(file.path(opt$ag_dir, sprintf("%s_GM12878.csv", vid)),
                   show_col_types = FALSE)
    ag %>% filter(gene_name == gene, track_strand %in% c("+", "-")) %>%
      mutate(track_label = label_track(track_name, track_strand)) %>%
      transmute(variant_id = vid, family = fam, gene_name = gene,
                track_label, track_strand, raw_score, quantile_score)
  }) %>% ungroup()
write_tsv(tabular, file.path(opt$out_dir, "supp_fig_AG_per_track_data.tsv"))
cat(sprintf("Wrote %s\n", file.path(opt$out_dir, "supp_fig_AG_per_track_data.tsv")))
