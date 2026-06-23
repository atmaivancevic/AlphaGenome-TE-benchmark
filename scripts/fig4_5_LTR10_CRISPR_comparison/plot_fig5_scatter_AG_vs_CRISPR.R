#!/usr/bin/env Rscript
#
# Fig 5 per-enhancer scatter: AG-predicted vs CRISPRi-observed effect size
# (Fig 3 G/H aesthetic, coloured by CRISPR significance). Per panel (one variant):
# X = CRISPRi log2FC (Ivancevic 2024 DESeq2), Y = AG log2FC (GeneMaskLFCScorer
# raw_score, stranded RNA only). Dots: red = sig down, blue = sig up, grey = n.s.
# Spearman/Pearson (signed) over all "Both" rows -> legend + sidecar TSV.
# Input: Supp Table 10 ("Both" rows = AG-scored genes within +/-500 kb that are
# also in the CRISPRi table).
#
# Example usage:
# Rscript scripts/fig4_5_LTR10_CRISPR_comparison/plot_fig5_scatter_AG_vs_CRISPR.R \
#   --variant LTR10.ATG12 \
#   --supp10  supptables/supp_table_10_LTR10_CRISPR_AG_predictions.tsv \
#   --output  figures/FIG5_FINAL/LTR10.ATG12_AG_vs_CRISPR_scatter.pdf

suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(readr); library(ggplot2)
  library(ggrepel); library(tibble); library(tidyr)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--variant", default = NULL, help = "Variant ID (e.g. LTR10.ATG12)"),
  make_option("--supp10",  default = "supptables/supp_table_10_LTR10_CRISPR_AG_predictions.tsv"),
  make_option("--output",  default = NULL,
              help = "Output PDF path (defaults to figures/FIG5_FINAL/<variant>_AG_vs_CRISPR_scatter.pdf)"),
  make_option("--xlim",    default = NULL, help = "Comma-separated 'lo,hi' for x-axis"),
  make_option("--ylim",    default = NULL, help = "Comma-separated 'lo,hi' for y-axis")
)))
if (is.null(opt$variant)) stop("--variant is required (e.g. LTR10.ATG12)")
if (is.null(opt$output))  opt$output <- sprintf("figures/FIG5_FINAL/%s_AG_vs_CRISPR_scatter.pdf", opt$variant)
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

# Read supp 10 and forward-fill the per-variant identity column so we can
# subset to one variant block (the table has variant_id blank on subsequent
# rows by design).
s10 <- read_tsv(opt$supp10, show_col_types = FALSE,
                col_types = cols(.default = "c"))
s10 <- s10 %>% mutate(variant_id_ff = na_if(variant_id, "")) %>%
  fill(variant_id_ff, .direction = "down")
df <- s10 %>% filter(variant_id_ff == opt$variant)
if (nrow(df) == 0) stop(sprintf("No rows for --variant %s in %s", opt$variant, opt$supp10))

# Keep only "Both" rows (genes with both AG raw_score and CRISPR log2FC)
df <- df %>%
  mutate(AG_RNA_raw_score      = suppressWarnings(as.numeric(AG_RNA_raw_score)),
         crispr_log2FoldChange = suppressWarnings(as.numeric(crispr_log2FoldChange)),
         crispr_padj           = suppressWarnings(as.numeric(crispr_padj))) %>%
  filter(!is.na(AG_RNA_raw_score) & !is.na(crispr_log2FoldChange))

# Significance / direction categorisation
ALPHA <- 0.05
df <- df %>%
  mutate(sig = !is.na(crispr_padj) & crispr_padj < ALPHA,
         cat = case_when(
           sig & crispr_log2FoldChange < 0 ~ "down",
           sig & crispr_log2FoldChange > 0 ~ "up",
           TRUE                            ~ "ns"))

cat_palette <- c(down = "#FF0000", up = "#0000FF", ns = "#808080")
cat_size    <- c(down = 2.8,        up = 2.8,       ns = 1.8)
cat_alpha   <- c(down = 0.9,        up = 0.9,       ns = 0.85)

# Stats: Spearman + Pearson (signed) over ALL Both rows (Avsec 2026: no filter on
# the dependent variable). cor.test wrapped in tryCatch so n<3 panels report NA.
safe_cor <- function(x, y, method) {
  if (length(x) < 3) return(list(estimate = NA_real_, p.value = NA_real_))
  res <- suppressWarnings(tryCatch(cor.test(x, y, method = method),
                                   error = function(e) NULL))
  if (is.null(res)) return(list(estimate = NA_real_, p.value = NA_real_))
  list(estimate = unname(res$estimate), p.value = res$p.value)
}
sp <- safe_cor(df$crispr_log2FoldChange, df$AG_RNA_raw_score, "spearman")
pe <- safe_cor(df$crispr_log2FoldChange, df$AG_RNA_raw_score, "pearson")
stats <- list(
  variant_id          = opt$variant,
  n_total             = nrow(df),
  n_sig_down          = sum(df$cat == "down"),
  n_sig_up            = sum(df$cat == "up"),
  n_ns                = sum(df$cat == "ns"),
  spearman_rho_signed = sp$estimate,
  spearman_p          = sp$p.value,
  pearson_r_signed    = pe$estimate,
  pearson_p           = pe$p.value
)

# Stat labels (Fig 3 convention): rho/r to 2 dp; "(ns)" when p>0.05 (exact p in TSV).
fmt_p <- function(p) if (!is.na(p) && p < ALPHA) sprintf("'p =' ~ '%.2g'", p) else "'(ns)'"
rho_lbl <- if (is.na(stats$spearman_rho_signed)) {
  "rho == 'n.d. (n<3)'"
} else {
  sprintf("rho == %0.2f ~ %s", stats$spearman_rho_signed, fmt_p(stats$spearman_p))
}
r_lbl <- if (is.na(stats$pearson_r_signed)) {
  "italic(r) == 'n.d. (n<3)'"
} else {
  sprintf("italic(r) == %0.2f ~ %s", stats$pearson_r_signed, fmt_p(stats$pearson_p))
}
n_lbl   <- sprintf("italic(n) == %d ~ '(' ~ %d ~ 'sig,' ~ %d ~ 'n.s.)'",
                   stats$n_total, stats$n_sig_down + stats$n_sig_up, stats$n_ns)

# Plot: significant dots drawn on top of n.s. so labels surface.
df <- df %>% arrange(factor(cat, levels = c("ns","down","up")))
p <- ggplot(df, aes(x = crispr_log2FoldChange, y = AG_RNA_raw_score, colour = cat)) +
  geom_hline(yintercept = 0, colour = "grey88", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey88", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey55", linewidth = 0.4) +
  geom_point(aes(size = cat, alpha = cat)) +
  # Label all dots (sig + n.s.); grey labels for n.s. match the dots.
  geom_text_repel(data = df %>% filter(cat != "ns"),
                  aes(label = gene_name),
                  colour = "black", size = 3.2, fontface = "italic",
                  segment.colour = "black", segment.alpha = 0.5, segment.size = 0.3,
                  box.padding = 0.5, point.padding = 0.4,
                  max.overlaps = Inf, seed = 1) +
  geom_text_repel(data = df %>% filter(cat == "ns"),
                  aes(label = gene_name),
                  colour = "grey45", size = 3.2, fontface = "italic",
                  segment.colour = "grey55", segment.alpha = 0.5, segment.size = 0.3,
                  # Push n.s. labels up/right to clear the sig cluster near the origin.
                  nudge_x = 0.6, nudge_y = -0.10,
                  box.padding = 0.5, point.padding = 0.4,
                  max.overlaps = Inf, seed = 1) +
  scale_colour_manual(values = cat_palette, guide = "none") +
  scale_size_manual(values = cat_size, guide = "none") +
  scale_alpha_manual(values = cat_alpha, guide = "none") +
  # ρ / r top-LEFT; n bottom-right (matches Fig 3 panels)
  annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.6, parse = TRUE,
           size = 3.8, label = rho_lbl) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.10, vjust = 3.4, parse = TRUE,
           size = 3.8, label = r_lbl) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.8, parse = TRUE,
           size = 3.4, label = n_lbl) +
  labs(x = expression("Observed effect (CRISPRi log"[2] * " fold change)"),
       y = "Predicted effect (AG raw score)",
       title = opt$variant) +
  theme_classic(base_size = 11) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
        axis.title = element_text(size = 11))

# Optional axis overrides; expand=FALSE pins edges to the given limits.
if (!is.null(opt$xlim) || !is.null(opt$ylim)) {
  xl <- if (!is.null(opt$xlim)) as.numeric(strsplit(opt$xlim, ",")[[1]]) else NULL
  yl <- if (!is.null(opt$ylim)) as.numeric(strsplit(opt$ylim, ",")[[1]]) else NULL
  p <- p + coord_cartesian(xlim = xl, ylim = yl, expand = FALSE)
  if (!is.null(xl)) {
    x_breaks <- seq(floor(xl[1]), ceiling(xl[2]), by = 1)
    p <- p + scale_x_continuous(breaks = x_breaks)
  }
  if (!is.null(yl)) {
    # 0.1 steps for y (AG raw_score range is typically <1)
    y_breaks <- seq(floor(yl[1] * 10) / 10, ceiling(yl[2] * 10) / 10, by = 0.1)
    p <- p + scale_y_continuous(breaks = y_breaks)
  }
}

ggsave(opt$output, p, width = 5, height = 4.5, dpi = 300)
message(sprintf("Wrote %s", opt$output))

# Stats sidecar
stats_path <- sub("\\.pdf$", "_stats.tsv", opt$output)
write_tsv(tibble::enframe(stats, "metric", "value") %>%
            mutate(value = as.character(value)), stats_path)
message(sprintf("Wrote %s", stats_path))

# Brief stdout summary
cat(sprintf("\n%s (n=%d: %d sig-down, %d sig-up, %d n.s.)\n",
            opt$variant, stats$n_total, stats$n_sig_down, stats$n_sig_up, stats$n_ns))
cat(sprintf("  Spearman rho (signed) = %+.3f  p = %.3g\n",
            stats$spearman_rho_signed, stats$spearman_p))
cat(sprintf("  Pearson  r   (signed) = %+.3f  p = %.3g\n",
            stats$pearson_r_signed,   stats$pearson_p))
