#!/usr/bin/env Rscript
#
# Per-candidate eQTL boxplot renderer for Fig 3. One PDF per (variant, gene):
# INT(residual expression) by genotype dose (0/1/2) with jittered points and a
# dosage-mean trend line, titled with the eQTL beta/p/q and the AG raw/quantile
# score (from Supp Table 5; AG line omitted when there's no AG score). Reads the
# cached cohort inputs from eqtl_matrixeqtl_pipeline.R.
#
# Example usage:
# Rscript scripts/fig3_polymorphic_TE_eQTLs/plot_fig3_eqtl_boxplots.R --cohort MAGE
# Rscript scripts/fig3_polymorphic_TE_eQTLs/plot_fig3_eqtl_boxplots.R --cohort GEUVADIS

suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(readr); library(ggplot2); library(cowplot); library(ggtext)
})

# Concordance banner lookup: maps the AG_RNA_concordance_with_MAGE value to
# display text + fill colour (a coloured strip above the plot title).
banner_for <- function(concordance) {
  switch(concordance,
    "match_high_confidence"    = list(text = "HIGH-CONFIDENCE CONCORDANCE",  fill = "#14512A"),
    "match"                    = list(text = "CONCORDANT",                   fill = "#2A7A41"),
    "mismatch_high_confidence" = list(text = "HIGH-CONFIDENCE DISCORDANCE",  fill = "#8B1A1F"),
    "mismatch"                 = list(text = "LOW-CONFIDENCE DISCORDANCE",   fill = "#C04848"),
    "AG_RNA_not_available"     = list(text = "AG SCORE UNAVAILABLE",         fill = "#888888"),
    "MAGE_beta_not_available"  = list(text = "NOT TESTED IN MAGE",           fill = "#888888"),
    list(text = "—", fill = "#BBBBBB"))
}

opt <- parse_args(OptionParser(option_list = list(
  make_option("--cohort",        default = NULL, help = "MAGE or GEUVADIS"),
  make_option("--variant_class", default = "INS",
              help = "INS (default; reads Supp 5 + INS eQTL outputs) or DEL (reads Supp 8 + DEL eQTL outputs). Switches all I/O paths + boxplot title text."),
  make_option("--candidates",    default = NULL,
              help = "Candidate-picks TSV; default depends on --variant_class"),
  make_option("--supp5",         default = NULL,
              help = "Supp 5 (INS) or Supp 8 (DEL) — for AG_RNA_raw_score / AG_RNA_quantile_score lookup per (variant, gene). Default depends on --variant_class."),
  make_option("--out_dir",       default = NULL, help = "Output dir; default depends on --variant_class"),
  make_option("--only_variant",  default = NULL,
              help = "If set, render only the candidate whose variant_id matches"),
  make_option("--width",         default = 4.0,  help = "Figure width in inches"),
  make_option("--height",        default = 4.7,  help = "Figure height in inches")
)))
if (is.null(opt$cohort)) stop("--cohort {MAGE,GEUVADIS} required")
opt$cohort <- toupper(opt$cohort)
stopifnot(opt$cohort %in% c("MAGE","GEUVADIS"))
opt$variant_class <- toupper(opt$variant_class)
stopifnot(opt$variant_class %in% c("INS", "DEL"))

# Variant-class handling — DEL uses parallel results dirs + Supp 8 + a DEL output
# dir; at render time the only difference is the title verb (insertion/deletion).
suffix <- if (opt$variant_class == "INS") "" else "_DEL"
cohort_dir <- paste0(if (opt$cohort == "MAGE") "results/eqtl_matrixeqtl_MAGE260"
                                          else "results/eqtl_matrixeqtl_GEUVADIS121",
                     suffix)
if (is.null(opt$candidates)) opt$candidates <- file.path(cohort_dir, "candidate_picks_working_set.tsv")
if (is.null(opt$supp5)) opt$supp5 <- if (opt$variant_class == "INS") {
    "supptables/supp_table_5_all_variant_gene_pairs.tsv"
  } else {
    "supptables/supp_table_8_all_variant_gene_pairs_DEL.tsv"
  }
int_rds   <- file.path(cohort_dir, sprintf("expression_INT_residuals_%s.rds", opt$cohort))
dose_rds  <- file.path(cohort_dir, sprintf("genotype_dosage_%s.rds",          opt$cohort))
eqtl_tsv  <- file.path(cohort_dir, "all_polymorphicTE_eqtls.tsv")

if (is.null(opt$out_dir)) {
  out_dir <- file.path(paste0("figures/fig3_eqtl_matrixeqtl", suffix),
                       if (opt$cohort == "MAGE") "MAGE260" else "GEUVADIS121")
} else {
  out_dir <- opt$out_dir
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("=== Rendering Fig 3 boxplots: %s / %s ===\n", opt$variant_class, opt$cohort))
cat(sprintf("Inputs:  %s, %s, %s\n", int_rds, dose_rds, eqtl_tsv))
cat(sprintf("Supp:    %s\n", opt$supp5))
cat(sprintf("Output:  %s/\n\n", out_dir))

# ---- Load cached pipeline outputs ----
int_obj  <- readRDS(int_rds)
dose_obj <- readRDS(dose_rds)
int_mat   <- int_obj$int_mat                    # genes × donors
gene_v    <- rownames(int_mat)                  # versioned ENSG IDs
gene_uv   <- sub("\\.[0-9]+$", "", gene_v)      # strip version suffix
dose_mat  <- dose_obj$dose_mat                  # variants × donors
testable  <- dose_obj$testable                  # variant metadata
stopifnot(identical(colnames(int_mat), colnames(dose_mat)))

eqtl_df <- read_tsv(eqtl_tsv, show_col_types = FALSE)

# ---- Load Supp 5 for per-(variant, gene) AG_RNA_raw_score / AG_RNA_quantile_score ----
supp5 <- read_tsv(opt$supp5, show_col_types = FALSE) %>%
  select(variant_id, pos, gene_symbol, AG_RNA_raw_score, AG_RNA_quantile_score,
         AG_RNA_concordance_with_MAGE)
cat(sprintf("Loaded Supp 5: %d (variant, gene) pairs (%d with AG scores)\n",
            nrow(supp5), sum(!is.na(supp5$AG_RNA_raw_score))))

# Per-variant target-gene metadata for the strand-aware distance label (in gene
# body / intergenic; N kb up/downstream of TSS). GENCODE v46 coords.
gene_info_locked <- tibble::tribble(
  ~variant_id,        ~gene_strand, ~gene_start, ~gene_end,
  # INS exemplars (Fig 3 panels A/B/C + Fig 2 LTR5_Hs)
  "SvimAsm00107100",  "+",          43680679,    43856617,    # HSD17B12 (Alu)
  "SvimAsm00022857",  "-",          190504337,   190534722,   # NEMP2 (L1)
  "SvimAsm00107233",  "+",          47248299,    47269033,    # NR1H3 (SVA)
  "SvimAsm00060017",  "+",          32741390,    32747198,    # HLA-DQA2 (LTR5_Hs)
  # DEL exemplars (Fig 3 panels D/E/F) — GENCODE v46 coords
  "SvimAsm00083149",  "+",          54135240,    54148514,    # MRPL15 (Alu)
  "SvimAsm00028490",  "-",          41246598,    41962130,    # ULK4 (L1)
  "SvimAsm00061340",  "+",          56955106,    57027346,    # BEND6 (SVA)
)

# ---- Load candidate pick set ----
candidates <- read_tsv(opt$candidates, show_col_types = FALSE) %>%
  select(variant_id, target_gene, role)
if (!is.null(opt$only_variant)) {
  candidates <- candidates %>% filter(variant_id == opt$only_variant)
  if (nrow(candidates) == 0) {
    stop(sprintf("--only_variant %s not in candidates file", opt$only_variant))
  }
}

# Typography — Fig 3 paper-figure defaults (150 x 200 pt panels: 11 pt axis
# titles, 8 pt ticks).
PLOT_FONT   <- "Helvetica"
FS_TITLE    <- 10
FS_SUBTITLE <- 7
FS_AXIS     <- 9
FS_TICK     <- 7
FS_BANNER   <- 2.8
DOT_SIZE    <- 0.7
SMOOTH_LW   <- 0.5

n_made <- 0L; n_skipped <- 0L; rows_for_log <- list()
for (i in seq_len(nrow(candidates))) {
  vid     <- candidates$variant_id[i]
  gn      <- candidates$target_gene[i]
  role    <- candidates$role[i]

  # eQTL stats for this (variant, gene)
  hit <- eqtl_df %>% filter(variant_id == vid, gene_symbol == gn)
  if (nrow(hit) == 0) {
    cat(sprintf("[%d/%d] SKIP %s/%s — not in eQTL output (not testable in %s)\n",
                i, nrow(candidates), vid, gn, opt$cohort))
    n_skipped <- n_skipped + 1L
    next
  }
  hit <- hit[1,]
  fam     <- hit$family
  # SVLEN is positive for INS, negative for DEL — display absolute value
  ins_len <- abs(hit$SVLEN)
  beta    <- hit$beta
  pval    <- hit$p
  qval    <- hit$q
  ensg_v_target <- hit$gene_id
  # TE-gene distance — eQTL pipeline writes it as `distance` (variant_pos - gene_TSS)
  distance_bp <- hit$distance
  distance_kb <- distance_bp / 1000

  # AG scores for this (variant, gene) pair
  ag_hit <- supp5 %>% filter(variant_id == vid, gene_symbol == gn)
  ag_raw <- if (nrow(ag_hit) >= 1) ag_hit$AG_RNA_raw_score[1]      else NA_real_
  ag_qtl <- if (nrow(ag_hit) >= 1) ag_hit$AG_RNA_quantile_score[1] else NA_real_
  ag_concordance <- if (nrow(ag_hit) >= 1) ag_hit$AG_RNA_concordance_with_MAGE[1] else NA_character_

  # Get expression vector for the target gene (use versioned ID for exact match)
  gene_idx <- match(ensg_v_target, gene_v)
  if (is.na(gene_idx)) {
    cat(sprintf("[%d/%d] SKIP %s/%s — target gene not in INT matrix\n",
                i, nrow(candidates), vid, gn))
    n_skipped <- n_skipped + 1L
    next
  }
  expr <- int_mat[gene_idx, ]

  # Genotype dosage
  dose <- dose_mat[vid, ]

  d <- data.frame(donor = colnames(int_mat), dose = factor(dose, levels = c(0L, 1L, 2L)),
                  expr = expr)
  n_per <- table(d$dose)
  geno_labels <- c("0/0", "0/1", "1/1")
  x_labels <- sprintf("%s\n(n=%d)", geno_labels, as.integer(n_per[c("0","1","2")]))

  # Compose title / subtitle
  cohort_label <- if (opt$cohort == "MAGE") sprintf("MAGE cohort n=%d",     ncol(int_mat))
                  else                       sprintf("GEUVADIS cohort n=%d", ncol(int_mat))

  # Strand-aware distance label (in gene body / intergenic), falling back to the
  # genome-coord signed distance when the variant isn't in gene_info_locked.
  variant_pos <- ag_hit$pos[1]
  gi <- gene_info_locked %>% filter(variant_id == vid)
  if (nrow(gi) == 1) {
    intragenic <- (variant_pos >= gi$gene_start) && (variant_pos <= gi$gene_end)
    # For - strand genes, flip the genome-coord distance to express it in
    # transcription-direction terms (positive = downstream of TSS in gene direction).
    strand_aware_kb <- if (gi$gene_strand == "-") -distance_kb else distance_kb
    direction <- if (strand_aware_kb > 0) "downstream" else "upstream"
    gene_loc <- if (intragenic) "In gene body" else "Intergenic"
    dist_tag <- sprintf("%s, %.1f kb %s of gene TSS",
                        gene_loc, abs(strand_aware_kb), direction)
  } else {
    dist_tag <- sprintf("%+.1f kb to TSS", distance_kb)
  }
  # Bold two-line title via explicit <br>: "{family} insertion ({size} bp):" /
  # "eQTL for gene {gene}". Verb (insertion/deletion) from sign(SVLEN).
  vtype_word <- if (hit$SVLEN < 0) "deletion" else "insertion"
  title_str <- sprintf("<b>%s %s (%d bp):<br>eQTL for gene %s</b>",
                       fam, vtype_word, ins_len, gn)

  # Subtitle (Avsec 2026 Fig 4d style): two prefix-aligned lines, headline numbers
  # (beta, AG effect) bolded via ggtext as eye anchors. "beta" written ASCII since
  # pdf() can't render Greek without cairo. Predicted line tinted with the banner colour.
  bn <- banner_for(if (is.na(ag_concordance)) "" else ag_concordance)
  beta_md   <- sprintf("<b>%+.3f</b>", beta)
  ag_raw_md <- if (is.na(ag_raw)) "NA" else sprintf("%+.3f", ag_raw)
  ag_qtl_md <- if (is.na(ag_qtl)) "NA" else sprintf("%+.3f", ag_qtl)
  # Predicted (AG) line first, Observed (MAGE) second; each headline value (6 pt)
  # followed by its parenthetical (5 pt) — AG quantile, MAGE FDR q. Predicted line
  # bold + banner-coloured. (Raw p dropped; q is in Supp 4/7.)
  line1 <- sprintf("<span style='color:%s'><b>Predicted effect size: %s</b> <span style='font-size:6pt'>(quantile = %s)</span></span>",
                   bn$fill, ag_raw_md, ag_qtl_md)
  line2 <- sprintf("Observed effect size: %s <span style='font-size:6pt'>(q = %.2e)</span>",
                   beta_md, qval)
  subtitle_str <- sprintf("%s<br>%s", line1, line2)

  # Compute group means for trend-line overlay
  means <- aggregate(expr ~ dose, data = d, FUN = mean)

  # Per-family 3-shade dose gradient (light 0/0 -> dark 1/1); boxes stay grey so
  # family identity reads from dot colour only (same hues carry into the Panel D
  # scatter). Alu uses grey (background family, 84% of variants); rare families keep saturated hues.
  dose_palettes <- list(
    Alu     = c("0" = "#CFCFCF", "1" = "#888888", "2" = "#4D4D4D"),
    L1      = c("0" = "#C9BFE3", "1" = "#7E6BBF", "2" = "#4A3982"),
    SVA     = c("0" = "#B8DDB3", "1" = "#5BA75A", "2" = "#2F6B3A"),
    LTR5_Hs = c("0" = "#9FC1D6", "1" = "#3A7CA5", "2" = "#1A4F73")
  )
  dose_palette <- if (fam %in% names(dose_palettes)) dose_palettes[[fam]]
                  else c("0" = "#CCCCCC", "1" = "#888888", "2" = "#444444")

  p <- ggplot(d, aes(x = dose, y = expr)) +
    geom_boxplot(outlier.shape = NA, width = 0.65, colour = "grey25",
                 fill = "grey90", alpha = 0.7) +
    geom_jitter(aes(colour = dose), width = 0.20, alpha = 0.65, size = DOT_SIZE) +
    scale_colour_manual(values = dose_palette, guide = "none") +
    # Linear regression line through individual points (not connecting
    # dose means) — matches the reference eQTL panel style.
    geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, se = FALSE,
                colour = "black", linewidth = SMOOTH_LW) +
    scale_x_discrete(labels = x_labels) +
    labs(title = title_str, subtitle = subtitle_str,
         x = "TE genotype",
         y = sprintf("INT(%s expression)", gn)) +
    theme_classic(base_size = 11, base_family = PLOT_FONT) +
    theme(plot.title    = element_markdown(size = FS_TITLE,    lineheight = 1.25,
                                           hjust = 0, family = PLOT_FONT),
          plot.subtitle = element_markdown(size = FS_SUBTITLE, colour = "grey25",
                                           hjust = 0, lineheight = 1.25,
                                           margin = margin(b = 10),
                                           family = PLOT_FONT),
          plot.title.position = "plot",
          axis.title    = element_text(size = FS_AXIS, family = PLOT_FONT),
          axis.text     = element_text(size = FS_TICK, family = PLOT_FONT),
          panel.grid = element_blank())

  # Concordance banner: coloured strip above the title showing the AG-vs-MAGE call
  # from Supp 5 (auto-updates when candidates change).
  banner_p <- ggplot() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = bn$fill) +
    annotate("text", x = 0.5, y = 0.5, label = bn$text, colour = "white",
             fontface = "bold", size = FS_BANNER, family = PLOT_FONT) +
    theme_void() +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE)

  final <- cowplot::plot_grid(banner_p, p, ncol = 1, rel_heights = c(0.08, 1))

  out_pdf <- file.path(out_dir, sprintf("%s_%s.pdf", vid, gsub("[^A-Za-z0-9_-]", "_", gn)))
  ggsave(out_pdf, final, width = as.numeric(opt$width),
         height = as.numeric(opt$height), device = "pdf")
  cat(sprintf("[%d/%d] %s\n", i, nrow(candidates), basename(out_pdf)))
  n_made <- n_made + 1L
  rows_for_log[[length(rows_for_log)+1]] <- data.frame(
    variant_id = vid, gene = gn, family = fam, SVLEN = hit$SVLEN,
    distance_kb = distance_kb,
    beta = beta, p = pval, q = qval,
    AG_RNA_raw_score = ag_raw, AG_RNA_quantile_score = ag_qtl,
    role = role, file = basename(out_pdf))
}

# Save a small index of what was rendered (only if we made anything)
log_path <- file.path(out_dir, "rendered_index.tsv")
if (length(rows_for_log) > 0) {
  write_tsv(do.call(rbind, rows_for_log), log_path)
} else {
  cat("(no PDFs rendered, skipping index)\n")
}

cat(sprintf("\nDone. %d PDFs written, %d skipped. Index: %s\n",
            n_made, n_skipped, log_path))
