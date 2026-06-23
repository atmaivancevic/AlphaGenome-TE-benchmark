#!/usr/bin/env Rscript
#
# scripts/fig3_polymorphic_TE_eQTLs/plot_fig3_scatter_AG_vs_MAGE.R
#
# Fig 3 Panel D — AG-predicted vs MAGE-observed effect size (Avsec 2026
# Fig 4d precedent). Two modes:
#   --mode top_gene  → Supp 4 (1 row per variant; top MAGE cis-gene)
#   --mode all_pairs → Supp 5 (1 row per (variant, gene) pair in cis)
#
# Stats reported (4):
#   - Spearman rho (signed)   on (beta_MAGE, AG_RNA_raw_score)
#   - Spearman rho (unsigned) on (|beta_MAGE|, |AG_RNA_raw_score|)
#   - Pearson r   (signed)
#   - Pearson r   (unsigned)
# Signed = does AG predict effect direction? Unsigned = does AG predict
# effect magnitude?
#
# Outputs go to figures/FIG3_FINAL{,_DEL}/panel{D,H}_*.pdf / *.tsv (8-panel
# Fig 3: A-D for INS, E-H for DEL; D and H are the AG-vs-MAGE scatter panels)

suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(readr); library(ggplot2)
  library(ggrepel); library(tibble)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--mode", default = "top_gene",
              help = "top_gene (one row per variant) | all_pairs (one row per (variant, gene) pair)"),
  make_option("--variant_class", default = "INS",
              help = "INS (default; reads Supp 4/5, writes to FIG3_FINAL/SUPP_FOR_FIG3) or DEL (reads Supp 7/8, writes to FIG3_DEL/SUPP_FOR_FIG3_DEL)"),
  make_option("--supp4",      default = NULL,
              help = "Top-gene supp table; default depends on --variant_class (Supp 4 INS / Supp 7 DEL)"),
  make_option("--supp5",      default = NULL,
              help = "All-pairs supp table; default depends on --variant_class (Supp 5 INS / Supp 8 DEL)"),
  make_option("--width",      default = 5.5, help = "Figure width in inches"),
  make_option("--height",     default = 5.5, help = "Figure height in inches"),
  make_option("--output",     default = NULL,
              help = "If set, overrides the default out_pdf path"),
  make_option("--title",      default = NULL,
              help = "Optional plot title (e.g. 'Polymorphic TE insertions')"),
  make_option("--exemplars",  default = NULL,
              help = "Comma-separated variant_ids to highlight (default: all locked exemplars)")
)))
stopifnot(opt$mode %in% c("top_gene", "all_pairs"))
opt$variant_class <- toupper(opt$variant_class)
stopifnot(opt$variant_class %in% c("INS", "DEL"))

# Resolve defaults based on variant_class
if (is.null(opt$supp4)) opt$supp4 <- if (opt$variant_class == "INS") {
    "supptables/supp_table_4_top_gene_per_variant.tsv"
  } else {
    "supptables/supp_table_7_top_gene_per_variant_DEL.tsv"
  }
if (is.null(opt$supp5)) opt$supp5 <- if (opt$variant_class == "INS") {
    "supptables/supp_table_5_all_variant_gene_pairs.tsv"
  } else {
    "supptables/supp_table_8_all_variant_gene_pairs_DEL.tsv"
  }

# Family palette. Alu = grey because it's 84% of points and functions as
# visual "background"; saturated colours for the rare families (L1 purple,
# SVA green, LTR5_Hs blue) so they pop. The Alu colour matches the bulk
# Alu boxplot fill so the two panels are visually consistent.
fam_palette <- c(Alu = "#888888", L1 = "#7E6BBF", SVA = "#5BA75A", LTR5_Hs = "#3A7CA5")

# Locked exemplars to highlight + label.
# INS exemplars are the 3 Fig 3 A/B/C boxplot variants + the Fig 2 HLA-DQA2
# exemplar. DEL exemplars are the 3 Fig 3 D/E/F boxplot variants, picked
# 2026-05-12 from Panel H's q-significant outlier candidates.
locked_INS <- tribble(
  ~variant_id,        ~gene_symbol,  ~exemplar_label,
  "SvimAsm00107100",  "HSD17B12",    "Alu : HSD17B12",
  "SvimAsm00022857",  "NEMP2",       "L1 : NEMP2",
  "SvimAsm00107233",  "NR1H3",       "SVA : NR1H3",
  "SvimAsm00060017",  "HLA-DQA2",    "LTR5_Hs : HLA-DQA2"
)
locked_DEL <- tribble(
  ~variant_id,        ~gene_symbol,  ~exemplar_label,
  "SvimAsm00083149",  "MRPL15",      "Alu : MRPL15",
  "SvimAsm00028490",  "ULK4",        "L1 : ULK4",
  "SvimAsm00061340",  "BEND6",       "SVA : BEND6"
)
locked <- if (opt$variant_class == "INS") locked_INS else locked_DEL
if (!is.null(opt$exemplars)) {
  keep_ids <- trimws(strsplit(opt$exemplars, ",", fixed = TRUE)[[1]])
  locked <- locked %>% filter(variant_id %in% keep_ids)
}

# Mode-specific input + column naming. Panel letters per locked 8-panel
# Fig 3 layout (2026-05-12): A/B/C = INS boxplots; D/E/F = DEL boxplots;
# G = INS scatter; H = DEL scatter. The two scatters sit side-by-side at
# the bottom of the figure for direct visual comparison.
# INS and DEL outputs go to the SAME dirs (FIG3_FINAL, SUPP_FOR_FIG3);
# filenames disambiguate via panel letter.
vtype_word   <- if (opt$variant_class == "INS") "TE insertions" else "TE deletions"
panel_letter <- if (opt$variant_class == "INS") "A"             else "B"
if (opt$mode == "top_gene") {
  raw <- read_tsv(opt$supp4, show_col_types = FALSE)
  df_full <- raw %>% rename(gene_symbol = top_gene_symbol_MAGE)
  n_label_template <- sprintf("'n' == '%%s %s'", vtype_word)
  out_dir   <- "figures/FIG3_FINAL"
  out_pdf   <- file.path(out_dir, sprintf("panel%s_%s_AG_vs_MAGE_scatter.pdf", panel_letter, opt$variant_class))
  out_stats <- file.path(out_dir, sprintf("panel%s_%s_AG_vs_MAGE_stats.tsv",   panel_letter, opt$variant_class))
} else {
  df_full <- read_tsv(opt$supp5, show_col_types = FALSE)
  n_label_template <- "'n' == '%s pairs'"
  out_dir   <- "figures/SUPP_FOR_FIG3"
  out_pdf   <- file.path(out_dir, sprintf("supp%s_%s_AG_vs_MAGE_scatter_allpairs.pdf", panel_letter, opt$variant_class))
  out_stats <- file.path(out_dir, sprintf("supp%s_%s_AG_vs_MAGE_stats_allpairs.tsv",   panel_letter, opt$variant_class))
}

n_input <- nrow(df_full)
df <- df_full %>%
  filter(!is.na(beta_MAGE), !is.na(AG_RNA_raw_score)) %>%
  mutate(family = factor(family, levels = names(fam_palette)))
n_plotted <- nrow(df)

cat(sprintf("=== Mode: %s ===\n", opt$mode))
cat(sprintf("Input rows: %d\n", n_input))
cat(sprintf("After NA filter (beta_MAGE × AG_RNA_raw_score both non-NA): %d\n", n_plotted))
cat(sprintf("Dropped: %d\n", n_input - n_plotted))
cat("Per-family counts:\n"); print(table(df$family))

# Four headline correlations + p-values (vs H0: cor = 0). cor.test() with
# method="spearman" warns about ties; exact=FALSE switches to the asymptotic
# normal approximation (appropriate at n≫30 anyway) and silences the warning.
ctest <- function(x, y, method) {
  suppressWarnings(cor.test(x, y, method = method, exact = FALSE))
}
t_rho_s  <- ctest(df$beta_MAGE,      df$AG_RNA_raw_score,      "spearman")
t_rho_u  <- ctest(abs(df$beta_MAGE), abs(df$AG_RNA_raw_score), "spearman")
t_r_s    <- ctest(df$beta_MAGE,      df$AG_RNA_raw_score,      "pearson")
t_r_u    <- ctest(abs(df$beta_MAGE), abs(df$AG_RNA_raw_score), "pearson")
stats <- list(
  spearman_rho_signed   = unname(t_rho_s$estimate),
  spearman_rho_unsigned = unname(t_rho_u$estimate),
  pearson_r_signed      = unname(t_r_s$estimate),
  pearson_r_unsigned    = unname(t_r_u$estimate)
)
pvals <- list(
  spearman_rho_signed_p   = t_rho_s$p.value,
  spearman_rho_unsigned_p = t_rho_u$p.value,
  pearson_r_signed_p      = t_r_s$p.value,
  pearson_r_unsigned_p    = t_r_u$p.value
)
cat("\nHeadline stats (estimate, p vs H0=0):\n")
for (nm in names(stats)) {
  pv <- pvals[[paste0(nm, "_p")]]
  cat(sprintf("  %-22s = %+.3f  (p = %.2e)\n", nm, stats[[nm]], pv))
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stats_df <- tibble(
  metric = c(names(stats), names(pvals), "n_input", "n_plotted"),
  value  = c(unlist(stats), unlist(pvals), n_input, n_plotted)
)
write_tsv(stats_df, out_stats)
cat(sprintf("\nWrote stats: %s\n", out_stats))

# Exemplars (always match Supp 4's top-gene definitions). LTR5_Hs /
# HLA-DQA2 gets its label handled by a dedicated geom_text_repel layer
# below (no leader line, hugged close above the point); the other three
# use default repel behaviour.
# Empty `locked` is allowed (DEL exemplars not yet picked) — the scatter
# just renders without highlighted points + labels in that case.
if (nrow(locked) > 0) {
  hl <- df %>% inner_join(locked, by = c("variant_id", "gene_symbol"))
  cat(sprintf("Highlighted exemplars matched: %d / %d\n", nrow(hl), nrow(locked)))
} else {
  hl <- df %>% slice(0) %>% mutate(exemplar_label = character(0))
  cat("No locked exemplars defined for this variant_class — rendering scatter without highlights.\n")
}

# Z-order: Alu (grey background family) drawn at the bottom, then the
# rare families (L1 / SVA / LTR5_Hs) shuffled together on top. With
# Alu = 84% of points, a pure random shuffle would bury most rare-family
# points under grey; layering rare-on-top makes every L1/SVA/LTR5_Hs
# point visible while the dense Alu cloud is still clearly the bulk.
# Within each band the shuffle is fixed by seed for reproducibility.
set.seed(20260511)
df <- bind_rows(
  df %>% filter(family == "Alu")  %>% slice_sample(prop = 1),
  df %>% filter(family != "Alu") %>% slice_sample(prop = 1)
)

# Build legend labels with per-family n
fam_counts <- table(factor(df$family, levels = names(fam_palette)))
fam_legend_labels <- sprintf("%s (n=%d)", names(fam_palette),
                             as.integer(fam_counts[names(fam_palette)]))
names(fam_legend_labels) <- names(fam_palette)

# Stat-annotation labels (computed outside the ggplot+ chain so we don't
# accidentally hit R's parser interpretation of `+ ALPHA <- 0.05` as the
# nonexistent `+<-` operator). Keep ρ/r to 2 dp; replace p-value with
# "(ns)" when p > 0.05 to reduce visual clutter — full p in the stats TSV.
ALPHA <- 0.05
fmt_p <- function(p) if (p < ALPHA) sprintf("'p =' ~ '%.2g'", p) else "'(ns)'"
rho_lbl <- sprintf("rho == %0.2f ~ %s",
                   stats$spearman_rho_signed, fmt_p(pvals$spearman_rho_signed_p))
r_lbl   <- sprintf("italic(r) == %0.2f ~ %s",
                   stats$pearson_r_signed,    fmt_p(pvals$pearson_r_signed_p))

# Typography — Fig 3 paper-figure defaults.
# Tuned for 220 × 220 pt panels: 10 pt axis titles + 9 pt ticks readable,
# legend kept at 7 pt so it fits inside the plot without crowding the points.
PLOT_FONT    <- "Helvetica"
FS_BASE      <- 10
FS_TITLE     <- 11
FS_AXIS      <- 10
FS_TICK      <- 8
FS_LEG_TITLE <- 8
FS_LEG_TEXT  <- 6
FS_STATS     <- 3.2
FS_N         <- 3.2
FS_EXEMPLAR  <- 2.5
PT_BULK      <- 1.2
PT_HL        <- 3.2

# Plot ----
p <- ggplot(df, aes(x = beta_MAGE, y = AG_RNA_raw_score, colour = family)) +
  geom_hline(yintercept = 0, colour = "grey88", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey88", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey55", linewidth = 0.4) +
  geom_point(size = PT_BULK, alpha = 0.5) +
  # Exemplars: large family-coloured point (uses the same palette as bulk;
  # Alu/HSD17B12 reads grey but pops via size + bold label).
  geom_point(data = hl, aes(colour = family), size = PT_HL) +
  # Other three exemplars: default repel placement, no leader lines.
  geom_text_repel(data = hl %>% filter(variant_id != "SvimAsm00060017"),
                  aes(label = exemplar_label),
                  colour = "black", size = FS_EXEMPLAR, fontface = "bold",
                  family = PLOT_FONT,
                  segment.colour = NA, box.padding = 0.7, point.padding = 0.4,
                  seed = 1, max.overlaps = Inf) +
  # LTR5_Hs / HLA-DQA2 (the Fig 2 exemplar): label sits directly above the
  # dot, hugged close, no leader line.
  geom_text_repel(data = hl %>% filter(variant_id == "SvimAsm00060017"),
                  aes(label = exemplar_label),
                  colour = "black", size = FS_EXEMPLAR, fontface = "bold",
                  family = PLOT_FONT,
                  nudge_y = 0.045, nudge_x = 0,
                  segment.colour = NA, box.padding = 0.1, point.padding = 0.1,
                  direction = "x", seed = 1) +
  scale_colour_manual(values = fam_palette, name = "TE family",
                      labels = fam_legend_labels, drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.6, parse = TRUE,
           size = FS_STATS, family = PLOT_FONT, label = rho_lbl) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.10, vjust = 3.4, parse = TRUE,
           size = FS_STATS, family = PLOT_FONT, label = r_lbl) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.10, vjust = -0.8, parse = TRUE,
           size = FS_N, family = PLOT_FONT,
           label = sprintf(n_label_template, format(n_plotted, big.mark = ","))) +
  coord_cartesian(xlim = c(-1.2, 1.2), ylim = c(-0.6, 0.6), expand = TRUE) +
  labs(title = opt$title,
       x = expression("Observed effect size across LCLs (MAGE " * beta * ")"),
       y = "Predicted effect size in GM12878 (AG)") +
  theme_classic(base_size = FS_BASE, base_family = PLOT_FONT) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = FS_TITLE, family = PLOT_FONT,
                                  face = "bold", hjust = 0,
                                  margin = margin(b = 4)),
        legend.position = "right",
        legend.justification = "top",
        legend.background = element_rect(fill = NA, colour = NA),
        legend.key.size = unit(0.35, "cm"),
        legend.margin = margin(0, 0, 0, 2),
        legend.box.margin = margin(0, 0, 0, -4),
        legend.title = element_text(size = FS_LEG_TITLE, family = PLOT_FONT),
        legend.text  = element_text(size = FS_LEG_TEXT,  family = PLOT_FONT),
        axis.title   = element_text(size = FS_AXIS, family = PLOT_FONT),
        axis.text    = element_text(size = FS_TICK, family = PLOT_FONT))

if (!is.null(opt$output)) out_pdf <- opt$output
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
ggsave(out_pdf, p,
       width  = as.numeric(opt$width),
       height = as.numeric(opt$height),
       device = "pdf")
cat(sprintf("Wrote scatter: %s\n", out_pdf))
