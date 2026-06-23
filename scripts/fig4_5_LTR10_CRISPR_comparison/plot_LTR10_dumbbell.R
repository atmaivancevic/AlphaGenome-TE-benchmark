#!/usr/bin/env Rscript
#
# Paired-lollipop ("dumbbell") comparing CRISPRi-observed RNA log2FC vs AG-predicted
# RNA raw_score per measured gene at one LTR10 CRISPRi enhancer. Two side-by-side
# panels share a y-axis (genes ordered by CRISPR log2FC); separate x-scales since
# the two magnitudes aren't directly comparable. Input: Supp Table 10.
#
# Example usage:
# Rscript scripts/fig4_5_LTR10_CRISPR_comparison/plot_LTR10_dumbbell.R --variant LTR10.ATG12

suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(readr); library(ggplot2)
  library(tidyr); library(patchwork)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--supp10", default = "supptables/supp_table_10_LTR10_CRISPR_AG_predictions.tsv"),
  make_option("--variant", default = "LTR10.ATG12"),
  make_option("--output", default = NULL,
              help = "Output PDF path; auto-derived from --variant if NULL"),
  make_option("--always_include", default = "",
              help = "Comma-separated gene names to force-include even if neither CRISPR nor AG passes the significance threshold (e.g. namesake target genes like MCPH1 for LTR10.MCPH1)."),
  # CRISPR x-axis shared across the 6 enhancers; AG x-axis "auto" (per-enhancer
  # symmetric range) so the XRCC4 outlier doesn't crush the others.
  make_option("--crispr_xlim", default = "-3.5,0.5"),
  make_option("--ag_xlim", default = "auto",
              help = "AG x-limits as 'lo,hi' or 'auto' for per-enhancer symmetric"),
  make_option("--width", default = 9.29, type = "numeric"),
  make_option("--height", default = NA, type = "numeric",
              help = "Plot height (in); auto-scales by gene count if NA"),
  make_option("--compact", action = "store_true", default = FALSE,
              help = "Compact mode: drops title + panel-headers, smaller fonts, tighter point sizes for Fig 5 multi-panel layout"),
  make_option("--pad_gene_names", default = 0, type = "integer",
              help = "Right-pad gene names to this character width with leading spaces. Use this to equalise axis-text widths (and therefore panel widths) across multiple dumbbells with different longest-gene-name lengths.")
)))
if (is.null(opt$output)) {
  opt$output <- sprintf("figures/FIG5_FINAL/dumbbells/%s_dumbbell.pdf", opt$variant)
}
parse_lim <- function(s) as.numeric(strsplit(s, ",", fixed = TRUE)[[1]])
CRISPR_XLIM <- parse_lim(opt$crispr_xlim)
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

# Read supp 10 and forward-fill variant_id (table has variant_id on the
# first row of each block and blank on subsequent rows).
s10 <- read_tsv(opt$supp10, show_col_types = FALSE, col_types = cols(.default = "c"))
s10 <- s10 %>% mutate(variant_id_ff = na_if(variant_id, "")) %>%
  fill(variant_id_ff, .direction = "down")

# Rule A: gene is CRISPR-tested and significant in CRISPR (padj<0.05) or AG
# (|quantile|>0.9). Each panel colours by its own direction/significance.
ALPHA   <- 0.05
AG_QTHR <- 0.9
df <- s10 %>%
  filter(variant_id_ff == opt$variant) %>%
  mutate(crispr_log2FoldChange  = suppressWarnings(as.numeric(crispr_log2FoldChange)),
         crispr_padj            = suppressWarnings(as.numeric(crispr_padj)),
         AG_RNA_raw_score       = suppressWarnings(as.numeric(AG_RNA_raw_score)),
         AG_RNA_quantile_score  = suppressWarnings(as.numeric(AG_RNA_quantile_score))) %>%
  filter(!is.na(crispr_log2FoldChange)) %>%
  mutate(crispr_sig = !is.na(crispr_padj) & crispr_padj < ALPHA,
         ag_sig     = !is.na(AG_RNA_quantile_score) &
                      abs(AG_RNA_quantile_score) > AG_QTHR) %>%
  { keep_always <- trimws(strsplit(opt$always_include, ",", fixed = TRUE)[[1]])
    filter(., crispr_sig | ag_sig | (gene_name %in% keep_always)) } %>%
  mutate(
    cat_crispr = case_when(
      crispr_sig & crispr_log2FoldChange < 0 ~ "down",
      crispr_sig & crispr_log2FoldChange > 0 ~ "up",
      TRUE                                   ~ "ns"),
    cat_ag = case_when(
      is.na(AG_RNA_raw_score)             ~ "ns",
      ag_sig & AG_RNA_raw_score < 0       ~ "down",
      ag_sig & AG_RNA_raw_score > 0       ~ "up",
      TRUE                                ~ "ns"))

# Order genes by CRISPR log2FC (most negative at top); optional right-pad so panel
# widths match across dumbbells.
df <- df %>% arrange(crispr_log2FoldChange)
if (opt$pad_gene_names > 0) {
  df$gene_name <- formatC(as.character(df$gene_name),
                          width = opt$pad_gene_names, flag = "")
}
df <- df %>% mutate(gene_name = factor(gene_name, levels = gene_name))

# Direction-and-significance palette. Dark red for sig-down + dark blue
# for sig-up matches the locked assay/effect colour convention used
# across the paper figures; grey for ns/NA.
cat_palette <- c(down = "#8B1A1F", up = "#1F3A8B", ns = "#808080")

# AG x-limits: "auto" gives a per-enhancer symmetric range padded 15%.
if (opt$ag_xlim == "auto") {
  m <- max(abs(range(df$AG_RNA_raw_score, na.rm = TRUE, finite = TRUE)))
  if (!is.finite(m) || m == 0) m <- 0.1
  AG_XLIM <- c(-m * 1.15, m * 1.15)
} else {
  AG_XLIM <- parse_lim(opt$ag_xlim)
}

FONT <- "Helvetica"
# Compact mode scales fonts/points down for Fig 5 multi-panel layout.
if (opt$compact) {
  # Match the browser-shot compact mode: 8pt body, 7pt tick labels.
  # STATS_SIZE is in ggplot mm units (size * 2.83 ≈ pt); 2.83 → ~8pt.
  BASE_SIZE <- 8; TITLE_SIZE <- 8; AXIS_TITLE_SIZE <- 8
  AXIS_TEXT_SIZE <- 7; GENE_TEXT_SIZE <- 8; POINT_SIZE <- 2.0
  STATS_SIZE <- 2.12; PANEL_TITLE_SIZE <- 8   # ~6pt stats (size * 2.83 ≈ pt)
  FIG_TITLE_SIZE <- 0  # 0 = hide figure title
} else {
  BASE_SIZE <- 18; TITLE_SIZE <- 20; AXIS_TITLE_SIZE <- 18
  AXIS_TEXT_SIZE <- 16; GENE_TEXT_SIZE <- 16; POINT_SIZE <- 4.5
  STATS_SIZE <- 6; PANEL_TITLE_SIZE <- 20; FIG_TITLE_SIZE <- 22
}
base_theme <- theme_classic(base_size = BASE_SIZE, base_family = FONT) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size = PANEL_TITLE_SIZE,
                                  family = FONT),
        axis.title = element_text(size = AXIS_TITLE_SIZE, family = FONT),
        axis.text  = element_text(size = AXIS_TEXT_SIZE, family = FONT))

# Panel L: CRISPR-observed log2FC, coloured by CRISPR direction + sig.
pL <- ggplot(df, aes(x = crispr_log2FoldChange, y = gene_name, colour = cat_crispr)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = crispr_log2FoldChange,
                   y = gene_name, yend = gene_name), linewidth = 0.7) +
  geom_point(size = POINT_SIZE) +
  scale_colour_manual(values = cat_palette, guide = "none") +
  scale_x_continuous(limits = CRISPR_XLIM, breaks = scales::pretty_breaks(n = 5)) +
  labs(x = expression("CRISPRi log"[2] * " fold change"),
       y = NULL,
       title = "Experiment") +
  base_theme +
  theme(axis.text.y = element_text(face = "italic", size = GENE_TEXT_SIZE,
                                   family = FONT))

# Panel R: AG-predicted raw_score on the per-enhancer AG_XLIM. Rows
# with AG_RNA_raw_score NA (out-of-window genes like MAOB at KDM6A)
# leave the row blank — y-axis factor preserves the slot.
pR <- ggplot(df, aes(x = AG_RNA_raw_score, y = gene_name, colour = cat_ag)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = AG_RNA_raw_score,
                   y = gene_name, yend = gene_name),
               linewidth = 0.7, na.rm = TRUE) +
  geom_point(size = POINT_SIZE, na.rm = TRUE) +
  scale_colour_manual(values = cat_palette, guide = "none") +
  scale_x_continuous(limits = AG_XLIM, breaks = scales::pretty_breaks(n = 5)) +
  labs(x = "Predicted effect size",
       y = NULL,
       title = "Prediction") +
  base_theme +
  theme(axis.text.y  = element_blank(),
        axis.ticks.y = element_blank())

# Correlation across the paired (both non-NA) observations. If fewer
# than 3 paired rows the annotation is skipped entirely — no "NA"
# artifact on enhancers like KDM6A with only n_paired = 1.
df_paired <- df %>% filter(!is.na(crispr_log2FoldChange) & !is.na(AG_RNA_raw_score))
safe_cor <- function(x, y, method) {
  if (length(x) < 3) return(list(estimate = NA_real_, p.value = NA_real_))
  res <- suppressWarnings(tryCatch(cor.test(x, y, method = method),
                                   error = function(e) NULL))
  if (is.null(res)) return(list(estimate = NA_real_, p.value = NA_real_))
  list(estimate = unname(res$estimate), p.value = res$p.value)
}
sp <- safe_cor(df_paired$crispr_log2FoldChange, df_paired$AG_RNA_raw_score, "spearman")
pe <- safe_cor(df_paired$crispr_log2FoldChange, df_paired$AG_RNA_raw_score, "pearson")

if (!is.na(sp$estimate)) {
  # Plotmath renders `rho` as the Greek letter (the plain PDF device
  # chokes on Unicode ρ).
  fmt_p <- function(p) if (!is.na(p) && p < ALPHA) sprintf("'p =' ~ '%.2g'", p) else "'(ns)'"
  rho_lbl <- sprintf("rho == %0.2f ~ %s", sp$estimate, fmt_p(sp$p.value))
  r_lbl   <- sprintf("italic(r) == %0.2f ~ %s", pe$estimate, fmt_p(pe$p.value))
  pR <- pR +
    annotate("text", x = -Inf, y = Inf, hjust = -0.07, vjust = 1.6, parse = TRUE,
             size = STATS_SIZE, family = FONT, label = rho_lbl) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.10, vjust = 3.4, parse = TRUE,
             size = STATS_SIZE, family = FONT, label = r_lbl)
}

p <- pL + pR + plot_layout(widths = c(1, 1))
if (!opt$compact) {
  p <- p + plot_annotation(
    title = sprintf("CRISPR-validated target genes for %s enhancer", opt$variant),
    theme = theme(plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE,
                                            hjust = 0.5, family = FONT)))
}

# Auto-scale height by number of rows so 27-gene MEF2D doesn't squish and
# 2-gene FGF2/KDM6A doesn't stretch.
if (is.na(opt$height)) {
  opt$height <- max(3.5, 2.0 + 0.30 * nrow(df))
}
ggsave(opt$output, p, width = opt$width, height = opt$height, dpi = 300)
message(sprintf("Wrote %s", opt$output))
