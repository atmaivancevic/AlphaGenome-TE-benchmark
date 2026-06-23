#!/usr/bin/env Rscript
#
# scripts/generate_supp_tables/build_supp_tables_4_5_7_8_eqtl.R
#
# Builds two TSV tabs that together summarize the polymorphic-TE cis-eQTL
# results across MAGE-260 and GEUVADIS-121:
#
#   Supp Table 4: supp_table_4_top_gene_per_variant.tsv  (1,322 rows)
#     One row per polymorphic-TE variant in MAGE-260.
#     For each, reports the variant's top cis gene (smallest MAGE p)
#     and the matching GEUVADIS β/p/q for that same (variant, gene).
#     Variants in gene deserts (no cis genes within +/- 500 kb) appear
#     with NA top-gene fields.
#
#   Supp Table 5: supp_table_5_all_variant_gene_pairs.tsv  (~9,234 rows)
#     One row per (variant, gene) cis pair tested in MAGE-260.
#     For each, reports MAGE β/p/q AND GEUVADIS β/p/q.
#     Standard long-format eQTL summary stats table.
#
# Both tabs are sorted by TE family (LTR5_Hs -> SVA -> L1 -> Alu -> misc)
# and then by MAGE q ascending. Both include a "cross_cohort_concordance"
# column with values:
#   same_direction       — MAGE and GEU β have same sign
#   direction_flip       — MAGE and GEU β have opposite signs
#   variant_not_testable_in_GEU — variant fails GEU's >=5 hom + >=5 ref filter
#   gene_or_pair_absent_in_GEU — variant testable in GEU but the (variant, gene) pair
#                                is absent from GEU eQTL output (most likely cause: gene
#                                not in GEU's GENCODE v12 annotation; could also be
#                                expression-filter dropout)
#   gene_desert_no_cis_in_MAGE — variant in MAGE set but has no cis gene
#   effect_size_zero           — one cohort produced β = 0 exactly (rare)
#
# Usage:
#   Rscript scripts/generate_supp_tables/build_supp_tables_4_5_7_8_eqtl.R
#
# Outputs are written to supptables/.

suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(readr)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--variant_class", default = "INS",
              help = "INS (default; outputs Supp 4 + 5) or DEL (outputs Supp 7 + 8)"),
  make_option("--mage_eqtl",   default = NULL,
              help = "MAGE eQTL TSV; defaults based on --variant_class"),
  make_option("--geu_eqtl",    default = NULL,
              help = "GEUVADIS eQTL TSV; defaults based on --variant_class"),
  make_option("--mage_dosage", default = NULL,
              help = "MAGE genotype-dosage RDS; defaults based on --variant_class"),
  make_option("--geu_dosage",  default = NULL,
              help = "GEUVADIS genotype-dosage RDS; defaults based on --variant_class"),
  make_option("--out_dir",     default = "supptables")
)))
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# Variant-class plumbing — DEL uses parallel results dirs (results/..._DEL/)
# and produces Supp Tables 7/8 instead of 4/5. Both arms share a single
# `SVLEN` length column (positive for INS, negative for DEL) sourced
# verbatim from the genotyped BCF in the Supp 3/6 builders.
opt$variant_class <- toupper(opt$variant_class)
stopifnot(opt$variant_class %in% c("INS", "DEL"))
suffix    <- if (opt$variant_class == "INS") ""        else "_DEL"
out_4_name <- if (opt$variant_class == "INS") {
  "supp_table_4_top_gene_per_variant.tsv"
} else {
  "supp_table_7_top_gene_per_variant_DEL.tsv"
}
out_5_name <- if (opt$variant_class == "INS") {
  "supp_table_5_all_variant_gene_pairs.tsv"
} else {
  "supp_table_8_all_variant_gene_pairs_DEL.tsv"
}
if (is.null(opt$mage_eqtl))   opt$mage_eqtl   <- sprintf("results/eqtl_matrixeqtl_MAGE260%s/all_polymorphicTE_eqtls.tsv", suffix)
if (is.null(opt$geu_eqtl))    opt$geu_eqtl    <- sprintf("results/eqtl_matrixeqtl_GEUVADIS121%s/all_polymorphicTE_eqtls.tsv", suffix)
if (is.null(opt$mage_dosage)) opt$mage_dosage <- sprintf("results/eqtl_matrixeqtl_MAGE260%s/genotype_dosage_MAGE.rds", suffix)
if (is.null(opt$geu_dosage))  opt$geu_dosage  <- sprintf("results/eqtl_matrixeqtl_GEUVADIS121%s/genotype_dosage_GEUVADIS.rds", suffix)
cat(sprintf("Variant class: %s -> %s + %s\n", opt$variant_class, out_4_name, out_5_name))

# ---- Load inputs ----
mage <- read_tsv(opt$mage_eqtl, show_col_types = FALSE)
geu  <- read_tsv(opt$geu_eqtl,  show_col_types = FALSE)
mage_dose <- readRDS(opt$mage_dosage)
geu_dose  <- readRDS(opt$geu_dosage)

mage_testable <- mage_dose$testable           # data.frame with columns variant_id, family, chrom, pos, SVLEN, n_hom, n_het, n_ref, maf
geu_testable_ids <- geu_dose$testable$variant_id      # variants that passed GEU filter

cat(sprintf("Loaded:\n"))
cat(sprintf("  MAGE eQTL pairs:     %d\n", nrow(mage)))
cat(sprintf("  GEU eQTL pairs:      %d\n", nrow(geu)))
cat(sprintf("  MAGE testable variants: %d\n", nrow(mage_testable)))
cat(sprintf("  GEU testable variants:  %d\n", length(geu_testable_ids)))

# ---- Sanity checks -------------------------------
# (1) GEU (variant_id, gene_id_uv) must be unique; otherwise the left_join below
#     would silently inflate row counts.
geu_dup <- geu %>% count(variant_id, gene_id_uv) %>% filter(n > 1)
if (nrow(geu_dup) > 0) {
  cat(sprintf("ERROR: %d duplicate (variant_id, gene_id_uv) pairs in GEU output:\n", nrow(geu_dup)))
  print(head(geu_dup, 10))
  stop("Aborting: duplicate keys would inflate cross-cohort joins.")
}
# (2) Every variant in MAGE eQTL output must be in MAGE testable set.
stopifnot(all(unique(mage$variant_id) %in% mage_testable$variant_id))

# ---- Family ordering helper ----
family_order <- c("LTR5_Hs" = 1L, "SVA" = 2L, "L1" = 3L, "Alu" = 4L)
order_family <- function(fam) {
  ifelse(fam %in% names(family_order), family_order[fam], 5L)  # 5 = misc
}

# ---- Concordance classifier (used by both tabs) ----
# Pre-compute lookup of GEU pairs as (variant_id, gene_id_uv) -> row.
# Matching on UNVERSIONED ENSG (gene_id_uv) instead of gene_symbol avoids
# gene-symbol-collision artifacts. GENCODE v12 has cases (e.g., RGS5) where
# two distinct ENSG IDs share the same HGNC symbol; symbol-based joins
# would inflate row counts via duplicate matches.
geu_lookup <- geu %>% select(variant_id, gene_id_uv,
                             gene_ensg_GEU = gene_id,           # versioned ENSG from GEU GENCODE v12
                             beta_GEU = beta, p_GEU = p, q_GEU = q,
                             n_hom_GEU121, n_het_GEU121, n_ref_GEU121, maf_GEU121)

classify_concordance <- function(beta_MAGE, beta_GEU, variant_in_geu_testable) {
  if (is.na(beta_MAGE)) return("gene_desert_no_cis_in_MAGE")
  if (is.na(beta_GEU)) {
    if (!variant_in_geu_testable) return("variant_not_testable_in_GEU")
    return("gene_or_pair_absent_in_GEU")
  }
  if (beta_MAGE == 0 || beta_GEU == 0) return("effect_size_zero")
  if (sign(beta_MAGE) == sign(beta_GEU)) return("same_direction")
  return("direction_flip")
}

# ============================================================================
# TAB 4 — top gene per variant (1,322 rows)
# ============================================================================
cat("\n[4] Building top-gene-per-variant table...\n")

# Top hit per variant in MAGE
mage_top <- mage %>%
  group_by(variant_id) %>%
  arrange(p, gene_id, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(variant_id, family,
            chrom = gene_chrom, pos = NA_integer_,   # placeholder; fill from testable below
            SVLEN,
            top_gene_symbol_MAGE = gene_symbol,
            top_gene_ensg_MAGE   = gene_id,
            top_gene_ensg_uv     = gene_id_uv,         # for unambiguous GEU join
            top_gene_TSS         = gene_TSS,
            distance_kb_MAGE     = distance / 1000,
            n_cis_genes_MAGE     = NA_integer_,        # filled below
            n_sig_q05_MAGE       = NA_integer_,        # filled below
            beta_MAGE = beta, p_MAGE = p, q_MAGE = q,
            n_hom_MAGE260, n_het_MAGE260, n_ref_MAGE260, maf_MAGE260)

# Per-variant cis-gene counts and significance counts
n_cis_per_variant <- mage %>% count(variant_id, name = "n_cis_genes_MAGE")
n_sig_per_variant <- mage %>% filter(q < 0.05) %>%
  count(variant_id, name = "n_sig_q05_MAGE")
mage_top <- mage_top %>%
  select(-n_cis_genes_MAGE, -n_sig_q05_MAGE) %>%
  left_join(n_cis_per_variant, by = "variant_id") %>%
  left_join(n_sig_per_variant, by = "variant_id") %>%
  mutate(n_sig_q05_MAGE = ifelse(is.na(n_sig_q05_MAGE), 0L, n_sig_q05_MAGE))

# Fill chrom/pos from variant metadata (canonical, not gene-coord)
mage_top <- mage_top %>% select(-chrom, -pos) %>%
  left_join(mage_testable %>% select(variant_id, chrom, pos, strand), by = "variant_id")

# Add gene-desert variants (testable in MAGE but with no cis genes)
desert_variants <- setdiff(mage_testable$variant_id, mage_top$variant_id)
desert_rows <- mage_testable %>% filter(variant_id %in% desert_variants) %>%
  transmute(variant_id, family, chrom, pos, SVLEN, strand,
            top_gene_symbol_MAGE = NA_character_,
            top_gene_ensg_MAGE   = NA_character_,
            top_gene_ensg_uv     = NA_character_,
            top_gene_TSS         = NA_integer_,
            distance_kb_MAGE     = NA_real_,
            n_cis_genes_MAGE     = 0L,
            n_sig_q05_MAGE       = 0L,
            beta_MAGE = NA_real_, p_MAGE = NA_real_, q_MAGE = NA_real_,
            n_hom_MAGE260 = n_hom, n_het_MAGE260 = n_het,
            n_ref_MAGE260 = n_ref, maf_MAGE260 = maf)

tab4 <- bind_rows(mage_top, desert_rows)
cat(sprintf("  MAGE-side rows: %d (with cis) + %d (gene desert) = %d total\n",
            nrow(mage_top), nrow(desert_rows), nrow(tab4)))

# Join GEU stats on (variant_id, unversioned ENSG)
tab4 <- tab4 %>%
  left_join(geu_lookup, by = c("variant_id", "top_gene_ensg_uv" = "gene_id_uv")) %>%
  mutate(variant_in_geu_testable = variant_id %in% geu_testable_ids) %>%
  rowwise() %>%
  mutate(cross_cohort_concordance = classify_concordance(
            beta_MAGE, beta_GEU, variant_in_geu_testable)) %>%
  ungroup() %>%
  select(-variant_in_geu_testable)

# Sort by family group, then MAGE q ascending
tab4 <- tab4 %>%
  mutate(family_rank = order_family(family),
         q_for_sort  = ifelse(is.na(q_MAGE), Inf, q_MAGE)) %>%
  arrange(family_rank, q_for_sort) %>%
  select(-family_rank, -q_for_sort)

# Final column order
tab4 <- tab4 %>%
  rename(top_gene_ensg_GEU = gene_ensg_GEU) %>%
  # Explicit FALSE for all non-significant cases, including missing q-values
  # (gene deserts on MAGE side; not-testable / absent pairs on GEU side).
  # The reason for missingness is preserved in cross_cohort_concordance.
  # (gene-desert variants on MAGE side; not-testable / absent on GEU side).
  mutate(top_gene_significant_MAGE = !is.na(q_MAGE) & q_MAGE < 0.05,
         top_gene_significant_GEU  = !is.na(q_GEU)  & q_GEU  < 0.05) %>%
  select(
    variant_id, family, chrom, pos, SVLEN, strand,
    n_hom_MAGE260, n_het_MAGE260, n_ref_MAGE260, maf_MAGE260,
    n_cis_genes_MAGE, n_sig_q05_MAGE,
    top_gene_symbol_MAGE, top_gene_ensg_MAGE, top_gene_TSS,
    distance_kb_MAGE, beta_MAGE, p_MAGE, q_MAGE, top_gene_significant_MAGE,
    top_gene_ensg_GEU,
    n_hom_GEU121, n_het_GEU121, n_ref_GEU121, maf_GEU121,
    beta_GEU, p_GEU, q_GEU, top_gene_significant_GEU,
    cross_cohort_concordance
  )
# Note: top_gene_ensg_uv was used as the join key but dropped from final output.
# Both top_gene_ensg_MAGE (v38) and top_gene_ensg_GEU (v12) are kept for traceability
# back to each cohort's native expression matrix.

out_4 <- file.path(opt$out_dir, out_4_name)
write_tsv(tab4, out_4)
cat(sprintf("  Saved %s (%d rows)\n", out_4, nrow(tab4)))

# ============================================================================
# TAB 5 — all variant-gene pairs (long format, ~9,234 rows)
# ============================================================================
cat("\n[5] Building all-variant-gene-pairs table...\n")

tab5 <- mage %>%
  transmute(variant_id, family,
            chrom = NA_character_, pos = NA_integer_,
            SVLEN, strand = NA_character_,
            gene_symbol, gene_ensg = gene_id, gene_ensg_uv = gene_id_uv, gene_TSS,
            distance_kb_MAGE = distance / 1000,
            beta_MAGE = beta, p_MAGE = p, q_MAGE = q,
            n_hom_MAGE260, n_het_MAGE260, n_ref_MAGE260, maf_MAGE260)

# Variant-level chrom / pos / strand
tab5 <- tab5 %>% select(-chrom, -pos, -strand) %>%
  left_join(mage_testable %>% select(variant_id, chrom, pos, strand), by = "variant_id")

# Join GEU stats on (variant_id, unversioned ENSG)
tab5 <- tab5 %>%
  left_join(geu_lookup, by = c("variant_id", "gene_ensg_uv" = "gene_id_uv")) %>%
  mutate(variant_in_geu_testable = variant_id %in% geu_testable_ids) %>%
  rowwise() %>%
  mutate(cross_cohort_concordance = classify_concordance(
            beta_MAGE, beta_GEU, variant_in_geu_testable)) %>%
  ungroup() %>%
  select(-variant_in_geu_testable)

# Sort: family group, then variant_id, then p_MAGE ascending within variant
tab5 <- tab5 %>%
  mutate(family_rank = order_family(family)) %>%
  arrange(family_rank, variant_id, p_MAGE, gene_ensg) %>%
  select(-family_rank)

tab5 <- tab5 %>%
  mutate(significant_q05_MAGE = !is.na(q_MAGE) & q_MAGE < 0.05,
         significant_q05_GEU  = !is.na(q_GEU)  & q_GEU  < 0.05) %>%
  select(
    variant_id, family, chrom, pos, SVLEN, strand,
    n_hom_MAGE260, n_het_MAGE260, n_ref_MAGE260, maf_MAGE260,
    gene_symbol, gene_ensg_MAGE = gene_ensg, gene_TSS, distance_kb_MAGE,
    beta_MAGE, p_MAGE, q_MAGE, significant_q05_MAGE,
    gene_ensg_GEU,
    n_hom_GEU121, n_het_GEU121, n_ref_GEU121, maf_GEU121,
    beta_GEU, p_GEU, q_GEU, significant_q05_GEU,
    cross_cohort_concordance
  )

out_5 <- file.path(opt$out_dir, out_5_name)
write_tsv(tab5, out_5)
cat(sprintf("  Saved %s (%d rows)\n", out_5, nrow(tab5)))

# ---- Summary ----
cat("\n=== Concordance breakdown ===\n")
cat("Tab 4 (per variant):\n")
print(table(tab4$cross_cohort_concordance, useNA = "ifany"))
cat("\nTab 5 (per variant-gene pair):\n")
print(table(tab5$cross_cohort_concordance, useNA = "ifany"))

cat("\n=== Family breakdown (Tab 4) ===\n")
print(tab4 %>% count(family, cross_cohort_concordance) %>% arrange(family, cross_cohort_concordance))
