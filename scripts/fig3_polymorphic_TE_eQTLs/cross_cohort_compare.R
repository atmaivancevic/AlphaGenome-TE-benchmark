#!/usr/bin/env Rscript
#
# Companion to eqtl_matrixeqtl_pipeline.R. Joins the MAGE-260 and GEUVADIS-121
# eQTL tables on (variant_id, gene_symbol) for a candidate pick set, producing a
# side-by-side cross-cohort table. The candidates TSV needs columns:
# variant_id, target_gene, role.
#
# Example usage:
# Rscript scripts/fig3_polymorphic_TE_eQTLs/cross_cohort_compare.R \
#   --candidates results/eqtl_matrixeqtl_MAGE260/candidate_picks_working_set.tsv \
#   --mage_eqtl  results/eqtl_matrixeqtl_MAGE260/all_polymorphicTE_eqtls.tsv \
#   --geu_eqtl   results/eqtl_matrixeqtl_GEUVADIS121/all_polymorphicTE_eqtls.tsv \
#   --out        results/eqtl_matrixeqtl_GEUVADIS121/candidate_picks_MAGE_vs_GEU.tsv

suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(readr)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--candidates", default = "results/eqtl_matrixeqtl_MAGE260/candidate_picks_working_set.tsv"),
  make_option("--mage_eqtl",  default = "results/eqtl_matrixeqtl_MAGE260/all_polymorphicTE_eqtls.tsv"),
  make_option("--geu_eqtl",   default = "results/eqtl_matrixeqtl_GEUVADIS121/all_polymorphicTE_eqtls.tsv"),
  make_option("--out",        default = "results/eqtl_matrixeqtl_GEUVADIS121/candidate_picks_MAGE_vs_GEU.tsv")
)))

cand <- read_tsv(opt$candidates, show_col_types = FALSE) %>%
  select(any_of(c("variant_id", "target_gene", "role")))   # minimal subset; MAGE/GEU stats are pulled fresh
mage <- read_tsv(opt$mage_eqtl,  show_col_types = FALSE)
geu  <- read_tsv(opt$geu_eqtl,   show_col_types = FALSE)

joined <- cand %>%
  left_join(mage %>% select(variant_id, gene_symbol,
                            family, ins_len, distance,
                            n_hom_MAGE260, n_het_MAGE260, n_ref_MAGE260, maf_MAGE260,
                            beta_MAGE = beta, p_MAGE = p, q_MAGE = q),
            by = c("variant_id","target_gene"="gene_symbol")) %>%
  left_join(geu %>% select(variant_id, gene_symbol,
                           n_hom_GEU121, n_het_GEU121, n_ref_GEU121, maf_GEU121,
                           beta_GEU = beta, p_GEU = p, q_GEU = q),
            by = c("variant_id","target_gene"="gene_symbol"))

# Replication status flag
joined <- joined %>%
  mutate(geu_status = case_when(
    is.na(beta_GEU)                      ~ "not_testable_in_GEU",
    sign(beta_MAGE) == sign(beta_GEU)    ~ "same_direction",
    TRUE                                  ~ "direction_flip"
  ))

write_tsv(joined, opt$out)
cat(sprintf("Saved %s (%d candidates)\n", opt$out, nrow(joined)))

# Print summary
cat("\n=== Cross-cohort summary ===\n")
cat(sprintf("%-19s %-13s %-9s %9s %9s | %9s %9s | %s\n",
            "variant_id","target","family","beta_MAGE","q_MAGE","beta_GEU","q_GEU","status"))
for (i in seq_len(nrow(joined))) {
  r <- joined[i,]
  bg <- if (is.na(r$beta_GEU)) "        --" else sprintf("%+9.3f", r$beta_GEU)
  qg <- if (is.na(r$q_GEU))    "        --" else sprintf("%9.2e", r$q_GEU)
  cat(sprintf("%-19s %-13s %-9s %+9.3f %9.2e | %s %s | %s\n",
              r$variant_id, r$target_gene, r$family,
              r$beta_MAGE, r$q_MAGE, bg, qg, r$geu_status))
}
