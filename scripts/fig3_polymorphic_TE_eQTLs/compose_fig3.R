#!/usr/bin/env Rscript
#
# scripts/fig3_polymorphic_TE_eQTLs/compose_fig3.R
#
# Compose Fig 3 panels A–D onto a single US Letter portrait, vector-PDF
# page suitable for opening in Illustrator. Layout:
#
#   ┌──────────┬──────────┬──────────┐
#   │  A  Alu  │  B  L1   │  C  SVA  │   ← 3 eQTL boxplots
#   └──────────┴──────────┴──────────┘
#   ┌──────────────────────────────────┐
#   │           D  AG vs MAGE          │   ← scatter
#   └──────────────────────────────────┘
#
# Self-contained: re-derives each panel's ggplot object (rather than
# embedding the already-rendered PDFs) so the output is fully vector,
# editable in Illustrator (text stays text, paths stay paths).

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(ggrepel)
  library(ggtext); library(cowplot); library(tibble)
})

# ──────────────────────────────────────────────────────────────────────────
# Common palette + banner mapper (must match plot_fig3_eqtl_boxplots.R and
# plot_fig3_scatter_AG_vs_MAGE.R)
# ──────────────────────────────────────────────────────────────────────────

fam_palette <- c(Alu = "#888888", L1 = "#7E6BBF", SVA = "#5BA75A", LTR5_Hs = "#3A7CA5")
dose_palettes <- list(
  Alu     = c("0" = "#CFCFCF", "1" = "#888888", "2" = "#4D4D4D"),
  L1      = c("0" = "#C9BFE3", "1" = "#7E6BBF", "2" = "#4A3982"),
  SVA     = c("0" = "#B8DDB3", "1" = "#5BA75A", "2" = "#2F6B3A"),
  LTR5_Hs = c("0" = "#9FC1D6", "1" = "#3A7CA5", "2" = "#1A4F73")
)

banner_for <- function(concordance) {
  switch(as.character(concordance),
    "match_high_confidence"    = list(text = "HIGH-CONFIDENCE CONCORDANCE",  fill = "#1A7F37"),
    "match"                    = list(text = "CONCORDANT",                   fill = "#4FAE5F"),
    "mismatch_high_confidence" = list(text = "HIGH-CONFIDENCE DISCORDANCE",  fill = "#B42318"),
    "mismatch"                 = list(text = "LOW-CONFIDENCE DISCORDANCE",   fill = "#DC6963"),
    list(text = "—", fill = "#BBBBBB"))
}

# Per-variant gene info for the strand-aware "In gene body / Intergenic"
# title line. GENCODE v46 coordinates.
gene_info_locked <- tribble(
  ~variant_id,        ~gene_strand, ~gene_start, ~gene_end,
  "SvimAsm00107100",  "+",          43680679,    43856617,
  "SvimAsm00022857",  "-",          190504337,   190534722,
  "SvimAsm00107233",  "+",          47248299,    47269033,
  "SvimAsm00060017",  "+",          32741390,    32747198,
)

# ──────────────────────────────────────────────────────────────────────────
# Load cohort + supp data once
# ──────────────────────────────────────────────────────────────────────────

cohort_dir <- "results/eqtl_matrixeqtl_MAGE260"
int_obj    <- readRDS(file.path(cohort_dir, "expression_INT_residuals_MAGE.rds"))
dose_obj   <- readRDS(file.path(cohort_dir, "genotype_dosage_MAGE.rds"))
eqtl_df    <- read_tsv(file.path(cohort_dir, "all_polymorphicTE_eqtls.tsv"), show_col_types = FALSE)
int_mat    <- int_obj$int_mat
gene_v     <- rownames(int_mat)
dose_mat   <- dose_obj$dose_mat
n_cohort   <- ncol(int_mat)

supp5 <- read_tsv("supptables/supp_table_5_all_variant_gene_pairs.tsv", show_col_types = FALSE) %>%
  select(variant_id, pos, gene_symbol, AG_RNA_raw_score, AG_RNA_quantile_score,
         AG_RNA_concordance_with_MAGE)

# ──────────────────────────────────────────────────────────────────────────
# Boxplot panel builder
# ──────────────────────────────────────────────────────────────────────────

build_boxplot_panel <- function(vid, gn) {
  hit <- eqtl_df %>% filter(variant_id == vid, gene_symbol == gn) %>% slice(1)
  fam <- hit$family; ins_len <- abs(hit$SVLEN)  # SVLEN is signed; title uses |SVLEN|
  beta <- hit$beta; pval <- hit$p; qval <- hit$q
  distance_kb <- hit$distance / 1000
  ensg_v_target <- hit$gene_id

  ag_hit <- supp5 %>% filter(variant_id == vid, gene_symbol == gn) %>% slice(1)
  ag_raw <- ag_hit$AG_RNA_raw_score
  ag_qtl <- ag_hit$AG_RNA_quantile_score
  ag_concordance <- ag_hit$AG_RNA_concordance_with_MAGE

  expr <- int_mat[ensg_v_target, ]
  dose <- dose_mat[vid, ]
  d <- data.frame(donor = colnames(int_mat),
                  dose = factor(dose, levels = c(0L, 1L, 2L)),
                  expr = expr)
  n_per <- table(d$dose)
  x_labels <- sprintf("%s\n(n=%d)", c("0/0", "0/1", "1/1"),
                      as.integer(n_per[c("0", "1", "2")]))

  variant_pos <- ag_hit$pos
  gi <- gene_info_locked %>% filter(variant_id == vid) %>% slice(1)
  intragenic <- (variant_pos >= gi$gene_start) && (variant_pos <= gi$gene_end)
  strand_aware_kb <- if (gi$gene_strand == "-") -distance_kb else distance_kb
  direction <- if (strand_aware_kb > 0) "downstream" else "upstream"
  gene_loc <- if (intragenic) "In gene body" else "Intergenic"
  dist_tag <- sprintf("%s, %.1f kb %s of gene TSS",
                      gene_loc, abs(strand_aware_kb), direction)

  title_str <- sprintf(
    "<b>%s insertion (%d bp): eQTL for %s</b><br><span style='font-weight:normal'>%s</span>",
    fam, ins_len, gn, dist_tag)

  beta_md   <- sprintf("<b>%+.3f</b>", beta)
  ag_raw_md <- sprintf("<b>%+.3f</b>", ag_raw)
  ag_qtl_md <- sprintf("%+.3f", ag_qtl)
  line1 <- sprintf("Observed effect size: %s  (p = %.2e, q = %.2e)",
                   beta_md, pval, qval)
  line2 <- sprintf("Predicted effect size: %s  (quantile = %s)",
                   ag_raw_md, ag_qtl_md)
  subtitle_str <- sprintf("<span style='font-size:8pt'>%s<br>%s</span>", line1, line2)

  dose_palette <- dose_palettes[[fam]]

  p <- ggplot(d, aes(x = dose, y = expr)) +
    geom_boxplot(outlier.shape = NA, width = 0.65, colour = "grey25",
                 fill = "grey90", alpha = 0.7) +
    geom_jitter(aes(colour = dose), width = 0.20, alpha = 0.65, size = 1.6) +
    scale_colour_manual(values = dose_palette, guide = "none") +
    geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, se = FALSE,
                colour = "black", linewidth = 0.5) +
    scale_x_discrete(labels = x_labels) +
    labs(title = title_str, subtitle = subtitle_str,
         x = "Polymorphic TE genotype",
         y = sprintf("INT(%s expression)", gn)) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_markdown(size = 10.5, lineheight = 1.25),
          plot.subtitle = element_markdown(size = 8, colour = "grey25",
                                           lineheight = 1.25),
          axis.title = element_text(size = 10),
          panel.grid = element_blank())

  bn <- banner_for(ag_concordance)
  banner_p <- ggplot() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = bn$fill) +
    annotate("text", x = 0.5, y = 0.5, label = bn$text, colour = "white",
             fontface = "bold", size = 3.5) +
    theme_void() +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE)

  plot_grid(banner_p, p, ncol = 1, rel_heights = c(0.06, 1))
}

# ──────────────────────────────────────────────────────────────────────────
# Scatter panel builder
# ──────────────────────────────────────────────────────────────────────────

build_scatter_panel <- function() {
  supp4 <- read_tsv("supptables/supp_table_4_top_gene_per_variant.tsv", show_col_types = FALSE)
  df_full <- supp4 %>% rename(gene_symbol = top_gene_symbol_MAGE)

  df <- df_full %>%
    filter(!is.na(beta_MAGE), !is.na(AG_RNA_raw_score)) %>%
    mutate(family = factor(family, levels = names(fam_palette)))
  n_plotted <- nrow(df)

  locked <- tribble(
    ~variant_id,        ~gene_symbol,  ~exemplar_label,
    "SvimAsm00107100",  "HSD17B12",    "Alu : HSD17B12",
    "SvimAsm00022857",  "NEMP2",       "L1 : NEMP2",
    "SvimAsm00107233",  "NR1H3",       "SVA : NR1H3",
    "SvimAsm00060017",  "HLA-DQA2",    "LTR5_Hs : HLA-DQA2"
  )
  hl <- df %>% inner_join(locked, by = c("variant_id", "gene_symbol"))

  ctest <- function(x, y, method) {
    suppressWarnings(cor.test(x, y, method = method, exact = FALSE))
  }
  t_rho_s <- ctest(df$beta_MAGE, df$AG_RNA_raw_score, "spearman")
  t_r_s   <- ctest(df$beta_MAGE, df$AG_RNA_raw_score, "pearson")

  set.seed(20260511)
  df <- bind_rows(
    df %>% filter(family == "Alu")  %>% slice_sample(prop = 1),
    df %>% filter(family != "Alu") %>% slice_sample(prop = 1)
  )

  fam_counts <- table(factor(df$family, levels = names(fam_palette)))
  fam_legend_labels <- sprintf("%s (n=%d)", names(fam_palette),
                               as.integer(fam_counts[names(fam_palette)]))

  ggplot(df, aes(x = beta_MAGE, y = AG_RNA_raw_score, colour = family)) +
    geom_hline(yintercept = 0, colour = "grey88", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey88", linewidth = 0.3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "grey55", linewidth = 0.4) +
    geom_point(size = 1.4, alpha = 0.5) +
    geom_point(data = hl, aes(colour = family), size = 4.5) +
    geom_text_repel(data = hl %>% filter(variant_id != "SvimAsm00060017"),
                    aes(label = exemplar_label),
                    colour = "black", size = 3.2, fontface = "bold",
                    segment.colour = NA, box.padding = 0.7, point.padding = 0.4,
                    seed = 1, max.overlaps = Inf) +
    geom_text_repel(data = hl %>% filter(variant_id == "SvimAsm00060017"),
                    aes(label = exemplar_label),
                    colour = "black", size = 3.2, fontface = "bold",
                    nudge_y = 0.045, nudge_x = 0,
                    segment.colour = NA, box.padding = 0.1, point.padding = 0.1,
                    direction = "x", seed = 1) +
    scale_colour_manual(values = fam_palette, name = "TE family",
                        labels = fam_legend_labels, drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(size = 2.0, alpha = 1))) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.6, parse = TRUE,
             size = 3.8, label = sprintf("rho == %0.2f * ',' ~ italic(p) == '%.1e'",
                                         t_rho_s$estimate, t_rho_s$p.value)) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.10, vjust = 3.4, parse = TRUE,
             size = 3.8, label = sprintf("italic(r) == %0.2f * ',' ~ italic(p) == '%.1e'",
                                         t_r_s$estimate, t_r_s$p.value)) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.10, vjust = -0.8, parse = TRUE,
             size = 4.0, label = sprintf("italic(n) == '%s TE insertions'",
                                          format(n_plotted, big.mark = ","))) +
    labs(x = expression("Observed effect size (MAGE " * beta * ")"),
         y = "Predicted effect size (AG raw score)") +
    theme_classic(base_size = 11) +
    theme(panel.grid = element_blank(),
          legend.position = c(0.97, 0.97),
          legend.justification = c(1, 1),
          legend.background = element_rect(fill = alpha("white", 0.8), colour = NA),
          legend.key.size = unit(0.4, "cm"),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 8),
          axis.title = element_text(size = 11))
}

# ──────────────────────────────────────────────────────────────────────────
# Build + compose
# ──────────────────────────────────────────────────────────────────────────

cat("Building panel A (Alu / HSD17B12)...\n")
pA <- build_boxplot_panel("SvimAsm00107100", "HSD17B12")
cat("Building panel B (L1 / NEMP2)...\n")
pB <- build_boxplot_panel("SvimAsm00022857", "NEMP2")
cat("Building panel C (SVA / NR1H3)...\n")
pC <- build_boxplot_panel("SvimAsm00107233", "NR1H3")
cat("Building panel D (AG vs MAGE scatter)...\n")
pD <- build_scatter_panel()

cat("Composing...\n")
top_row <- plot_grid(pA, pB, pC, ncol = 3,
                     labels = c("A", "B", "C"),
                     label_size = 14, label_fontface = "bold")
fig3 <- plot_grid(top_row, pD, ncol = 1,
                  rel_heights = c(1.0, 1.4),
                  labels = c("", "D"),
                  label_size = 14, label_fontface = "bold")

out_path <- "figures/FIG3_FINAL/Fig3.pdf"
ggsave(out_path, fig3, width = 8.5, height = 11, units = "in",
       device = "pdf")
cat(sprintf("\nWrote: %s (US Letter portrait, vector PDF — Illustrator-editable)\n", out_path))
