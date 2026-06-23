#!/usr/bin/env Rscript
#
# scripts/fig3_polymorphic_TE_eQTLs/eqtl_matrixeqtl_pipeline.R
#
# Polymorphic-TE cis-eQTL pipeline for Fig 3 (Ivancevic et al.).
# Runs MatrixEQTL on either MAGE-260 or GEUVADIS-121 with a unified covariate
# framework and outputs the full per-(variant, gene) table for downstream
# plotting + cross-cohort comparison.
#
# Pipeline (matches Bravo et al. 2024 LCL polymorphic-TE eQTL methodology):
#   1. Load cohort-specific expression matrix (MAGE: precomputed VST;
#      GEUVADIS: precomputed RPKM, log2(RPKM+1), aggregate sample-runs to
#      per-donor mean).
#   2. Restrict to Taylor et al. 2024's expression-quality-filtered gene set
#      (unversioned ENSG match across GENCODE versions where needed).
#   3. Residualize per gene against cohort-specific covariates:
#        MAGE:     INT(VST) ~ continentalGroup + sex + batch
#        GEUVADIS: INT(log2(RPKM+1)) ~ continentalGroup + sex
#      Inverse-normal-transform the residuals per gene
#      (Beasley, Erickson & Allison 2009 formula).
#   4. Build polymorphic-TE genotype dosage matrix from Schloissnig 2025
#      carrier lists (Supp Table 3). Filter: NA12878 = 0/0 + MAF > 0.01 +
#      >= 5 hom-alt + >= 5 hom-ref carriers in the cohort.
#   5. Run MatrixEQTL: additive linear model on INT(residuals) ~ dose,
#      cisDist = +/- 500 kb (matches AlphaGenome's 1 Mb receptive field),
#      pvOutputThreshold = 1.0 (save all (variant, gene) pairs for the
#      genome-wide AG-vs-observed scatter; multiple-testing in R post-hoc).
#   6. BH-FDR across all (variant, gene) pairs.
#   7. Annotate with gene/variant metadata and save TSV.
#
# Citations chain:
#   Bravo et al. 2024 (PLOS Genetics; LCL polymorphic-TE eQTL methodology;
#   covariate residualization-then-INT order)
#   Castanera et al. 2023 (eLife; subpopulation as covariate, no PEER)
#   Koks et al. 2021 (IJMS; sparse-covariate TE-eQTL precedent)
#   Jia et al. 2025 (review; "standard eQTL identification = linear regression")
#   Shabalin 2012 (Bioinformatics; MatrixEQTL)
#   Beasley, Erickson & Allison 2009 (Behavior Genetics; INT formula)
#   Taylor et al. 2024 (Nature; MAGE expression matrix and gene-filter list)
#
# Usage:
#   Rscript scripts/fig3_polymorphic_TE_eQTLs/eqtl_matrixeqtl_pipeline.R --cohort MAGE
#   Rscript scripts/fig3_polymorphic_TE_eQTLs/eqtl_matrixeqtl_pipeline.R --cohort GEUVADIS
#
# Outputs go to results/eqtl_matrixeqtl_{MAGE260,GEUVADIS121}/ by default.

suppressPackageStartupMessages({
  library(optparse); library(data.table); library(R.utils)
  library(dplyr); library(readr); library(matrixStats); library(MatrixEQTL)
})

opt_list <- list(
  make_option("--cohort", type = "character", help = "MAGE or GEUVADIS"),
  make_option("--variant_class", default = "INS",
              help = "INS (default; uses Supp Table 3) or DEL (Supp Table 6)"),
  make_option("--supp_table_3",  default = "supptables/supp_table_3_polymorphic_TE_insertions.tsv",
              help = "INS catalog (used when --variant_class INS)"),
  make_option("--supp_table_6",  default = "supptables/supp_table_6_polymorphic_TE_deletions.tsv",
              help = "DEL catalog (used when --variant_class DEL)"),
  make_option("--supp_table_2",  default = "supptables/supp_table_2_individuals.tsv"),
  make_option("--mage_vst",      default = "data/rna_cohorts/precomputed/MAGE/global_trend_results/global_expression_trends/expression.vst.csv"),
  make_option("--mage_metadata", default = "data/rna_cohorts/precomputed/MAGE/sample_library_info/sample.metadata.MAGE.v1.0.txt"),
  make_option("--mage_genelist", default = "data/rna_cohorts/precomputed/MAGE/QTL_results/eQTL_results/expression_filteredGenes.MAGE.v1.0.txt.gz"),
  make_option("--mage_overlap",  default = "data/rna_cohorts/schloissnig_phased_x_mage.tsv"),
  make_option("--geu_rpkm",      default = "data/rna_cohorts/precomputed/GEUVADIS/GD660.GeneQuantRPKM.txt.gz"),
  make_option("--maf_floor",     default = 0.01,    type = "double"),
  make_option("--min_hom",       default = 5L,      type = "integer"),
  make_option("--min_ref",       default = 5L,      type = "integer"),
  make_option("--cis_dist",      default = 500000L, type = "integer"),
  make_option("--out_dir",       default = NULL,
              help = "Output dir; if NULL, defaults to results/eqtl_matrixeqtl_<cohort_size>[_DEL]")
)
opt <- parse_args(OptionParser(option_list = opt_list))
if (is.null(opt$cohort)) stop("--cohort {MAGE,GEUVADIS} is required")
opt$cohort <- toupper(opt$cohort)
stopifnot(opt$cohort %in% c("MAGE","GEUVADIS"))
opt$variant_class <- toupper(opt$variant_class)
stopifnot(opt$variant_class %in% c("INS", "DEL"))

# Variant-class plumbing: INS uses Supp Table 3; DEL uses Supp Table 6.
# Both catalogs use a single `SVLEN` column (positive for INS, negative for
# DEL — sourced verbatim from the genotyped BCF's SVLEN INFO field). NA12878
# = 0/0 filter is the same syntax in both cases (REF-unbiased AG prediction),
# with opposite biology — INS carriers gained the TE; DEL carriers lost it.
catalog_path   <- if (opt$variant_class == "INS") opt$supp_table_3 else opt$supp_table_6
out_dir_suffix <- if (opt$variant_class == "INS") ""               else "_DEL"

# Inverse normal transform per Beasley et al. 2009
INT <- function(x) qnorm((rank(x, ties.method = "average") - 0.5) / length(x))

# Defensive: MAGE-published sample-swap pair (HG00237 / NA11919). Neither is in
# our 260 overlap, but defence-in-depth in case the metadata file changes.
SAMPLE_SWAP_SKIP <- c("HG00237", "NA11919")

# ---- Taylor's expression-quality-filtered gene set (used by both cohorts) ----
taylor <- read_tsv(opt$mage_genelist, show_col_types = FALSE) %>%
  mutate(ensg_uv = sub("\\.[0-9]+$", "", ensemblID),
         tss     = ifelse(strand == "+", chromStart + 1L, chromEnd))
cat(sprintf("[%s] Taylor expression-filtered gene set: %d genes\n", opt$cohort, nrow(taylor)))

# ============================================================================
# COHORT-SPECIFIC EXPRESSION LOADING + DESIGN
# ============================================================================
if (opt$cohort == "MAGE") {

  cohort <- read_tsv(opt$mage_overlap, show_col_types = FALSE)
  donors_target <- setdiff(cohort$sample, SAMPLE_SWAP_SKIP)

  # .keep_all = TRUE retains the per-library continuous covariates
  # (RIN, numReads, RNA Qubit fields) that the PC diagnostic needs. Picks the
  # first library per donor — donor-level fields (sex, super-pop, etc.) are
  # invariant across libraries; continuous fields are library-level and
  # represent that one library's value.
  mage_md <- read.table(opt$mage_metadata, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE) %>%
    distinct(sample_coriellID, .keep_all = TRUE)

  cat(sprintf("[MAGE] Loading VST matrix...\n"))
  vst_dt <- fread(opt$mage_vst)
  ensg_v <- vst_dt[[1]]
  vst <- as.data.frame(vst_dt[, -1]); rownames(vst) <- ensg_v

  # Rename coriell IDs -> kgp_id and intersect with target donors
  coriell_to_kgp <- setNames(mage_md$sample_kgpID, mage_md$sample_coriellID)
  match_idx <- match(colnames(vst), names(coriell_to_kgp))
  vst <- vst[, !is.na(match_idx), drop = FALSE]
  colnames(vst) <- coriell_to_kgp[colnames(vst)]
  vst <- vst[, intersect(colnames(vst), donors_target), drop = FALSE]

  # Restrict to Taylor's filter list (versioned IDs match exactly within MAGE)
  keep_genes <- rownames(vst) %in% taylor$ensemblID
  expr_mat <- as.matrix(vst[keep_genes, , drop = FALSE])
  gene_versioned   <- rownames(expr_mat)
  gene_unversioned <- sub("\\.[0-9]+$", "", gene_versioned)

  donor_meta <- mage_md %>%
    rename(donor = sample_kgpID) %>%
    filter(donor %in% colnames(expr_mat)) %>%
    distinct(donor, .keep_all = TRUE) %>%
    arrange(match(donor, colnames(expr_mat))) %>%
    select(donor, sex, continentalGroup, batch)
  stopifnot(identical(donor_meta$donor, colnames(expr_mat)))

  donor_meta$continentalGroup <- factor(donor_meta$continentalGroup)
  donor_meta$sex              <- factor(donor_meta$sex)
  donor_meta$batch            <- factor(donor_meta$batch)
  cov_formula  <- ~ continentalGroup + sex + batch
  cohort_label <- sprintf("MAGE-%d", ncol(expr_mat))
  out_dir_default <- paste0("results/eqtl_matrixeqtl_MAGE260", out_dir_suffix)

} else {  # GEUVADIS

  st2 <- read_tsv(opt$supp_table_2, show_col_types = FALSE)
  geu_meta <- st2 %>% filter(tier1_rna_geuvadis) %>%
    select(donor = sample, superpop, sex) %>% arrange(donor)

  cat(sprintf("[GEU] Loading RPKM matrix...\n"))
  rpkm_raw <- fread(opt$geu_rpkm)
  sample_cols <- setdiff(colnames(rpkm_raw), c("TargetID","Gene_Symbol","Chr","Coord"))
  ensg_v <- rpkm_raw$TargetID
  ensg_uv <- sub("\\.[0-9]+$", "", ensg_v)

  # Aggregate sample-runs to per-donor mean (some donors have multiple sequencing
  # runs in GEUVADIS; standard Lappalainen 2013 preprocessing).
  run_to_donor <- sub("\\..*$", "", sample_cols)
  keep_donors  <- intersect(unique(run_to_donor), geu_meta$donor)
  rpkm_mat <- as.matrix(rpkm_raw[, ..sample_cols]); rownames(rpkm_mat) <- ensg_v
  agg <- matrix(NA_real_, nrow(rpkm_mat), length(keep_donors),
                dimnames = list(ensg_v, keep_donors))
  for (d in keep_donors) {
    cols <- which(run_to_donor == d)
    agg[, d] <- if (length(cols) == 1) rpkm_mat[, cols] else rowMeans(rpkm_mat[, cols])
  }
  log2_mat <- log2(agg + 1)

  # Cross-version GENCODE intersect: strip version suffix on both sides
  keep_genes <- ensg_uv %in% taylor$ensg_uv
  expr_mat <- log2_mat[keep_genes, , drop = FALSE]
  gene_versioned   <- ensg_v[keep_genes]
  gene_unversioned <- ensg_uv[keep_genes]

  donor_meta <- geu_meta %>% filter(donor %in% keep_donors) %>%
    arrange(match(donor, colnames(expr_mat))) %>%
    rename(continentalGroup = superpop)
  stopifnot(identical(donor_meta$donor, colnames(expr_mat)))

  donor_meta$continentalGroup <- factor(donor_meta$continentalGroup)
  donor_meta$sex              <- factor(donor_meta$sex)
  cov_formula  <- ~ continentalGroup + sex
  cohort_label <- sprintf("GEUVADIS-%d", ncol(expr_mat))
  out_dir_default <- paste0("results/eqtl_matrixeqtl_GEUVADIS121", out_dir_suffix)
}

if (is.null(opt$out_dir)) opt$out_dir <- out_dir_default
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("\n[%s] Cohort: %s, %d genes x %d donors\n",
            opt$cohort, cohort_label, nrow(expr_mat), ncol(expr_mat)))
cat("  superpop:"); print(table(donor_meta$continentalGroup))
cat("  sex:");      print(table(donor_meta$sex))
if ("batch" %in% colnames(donor_meta))
  cat(sprintf("  batch: %d levels\n", nlevels(donor_meta$batch)))

# ============================================================================
# PRE-FLIGHT EXPRESSION-PC x COVARIATE DIAGNOSTIC (MAGE only; -> Supp Table 9)
# Top 10 PCs of the expression-filtered VST matrix; for each PC, the variance
# explained (R^2 from a univariate linear model) by every candidate covariate.
# Justifies the sparse covariate set against alternatives (genotype PCs, PEER).
# ============================================================================
if (opt$cohort == "MAGE") {
  cat(sprintf("\n[MAGE] Computing PC x covariate diagnostic...\n"))
  pca <- prcomp(t(expr_mat), center = TRUE, scale. = FALSE)
  n_pc <- 10
  pc_pct_var <- (pca$sdev^2 / sum(pca$sdev^2)) * 100

  diag_md <- mage_md %>%
    rename(donor = sample_kgpID) %>%
    filter(donor %in% colnames(expr_mat)) %>%
    distinct(donor, .keep_all = TRUE) %>%
    arrange(match(donor, colnames(expr_mat)))
  stopifnot(identical(diag_md$donor, colnames(expr_mat)))

  cov_specs <- list(
    list(name = "batch",            type = "categorical", col = "batch"),
    list(name = "continentalGroup", type = "categorical", col = "continentalGroup"),
    list(name = "population",       type = "categorical", col = "population"),
    list(name = "sex",              type = "categorical", col = "sex"),
    list(name = "RIN",              type = "continuous",  col = "RIN"),
    list(name = "numReads",         type = "continuous",  col = "numReads"),
    list(name = "RNAQubitConc",     type = "continuous",  col = "RNAQubitConc_ng.ul"),
    list(name = "RNAQubitAmount",   type = "continuous",  col = "RNAQubitTotalAmount_ng")
  )

  diag_rows <- list()
  for (i in seq_len(n_pc)) {
    pc_vec <- pca$x[, i]
    for (cs in cov_specs) {
      x <- diag_md[[cs$col]]
      if (cs$type == "categorical") x <- factor(x)
      r2 <- summary(lm(pc_vec ~ x))$r.squared
      diag_rows[[length(diag_rows) + 1]] <- data.frame(
        PC = sprintf("PC%d", i), pct_var_PC = pc_pct_var[i],
        covariate = cs$name, type = cs$type, pct_variance = r2 * 100)
    }
  }
  diag_df <- do.call(rbind, diag_rows)

  diag_results <- file.path(opt$out_dir, "expression_PC_diagnostic.tsv")
  diag_supp    <- "supptables/supp_table_9_pc_covariate_diagnostic.tsv"
  dir.create(dirname(diag_supp), recursive = TRUE, showWarnings = FALSE)
  write_tsv(diag_df, diag_results)
  write_tsv(diag_df, diag_supp)
  cat(sprintf("[MAGE] PC diagnostic: %d rows -> %s + %s\n",
              nrow(diag_df), diag_results, diag_supp))
}

# ============================================================================
# RESIDUALIZE + INVERSE NORMAL TRANSFORM
# ============================================================================
X   <- model.matrix(cov_formula, data = donor_meta)
qrX <- qr(X)
cat(sprintf("\n[%s] Design matrix: %d x %d  rank=%d\n",
            opt$cohort, nrow(X), ncol(X), qrX$rank))
if (qrX$rank < ncol(X)) {
  collinear_cols <- setdiff(seq_len(ncol(X)), qrX$pivot[seq_len(qrX$rank)])
  stop("Design matrix is rank-deficient; collinear columns: ",
       paste(colnames(X)[collinear_cols], collapse = ", "))
}

resid_mat <- t(qr.resid(qrX, t(expr_mat)))
dimnames(resid_mat) <- dimnames(expr_mat)
stopifnot(!anyNA(resid_mat) && all(is.finite(resid_mat)))

zero_var <- which(rowVars(resid_mat) < 1e-10)
if (length(zero_var) > 0) {
  resid_mat        <- resid_mat[-zero_var, , drop = FALSE]
  gene_versioned   <- gene_versioned[-zero_var]
  gene_unversioned <- gene_unversioned[-zero_var]
}

int_mat <- t(apply(resid_mat, 1, INT))
dimnames(int_mat) <- dimnames(resid_mat)
cat(sprintf("[%s] INT matrix: %d genes x %d donors  (per-gene sd median=%.4f)\n",
            opt$cohort, nrow(int_mat), ncol(int_mat), median(rowSds(int_mat))))

saveRDS(list(int_mat = int_mat,
             gene_versioned = gene_versioned, gene_unversioned = gene_unversioned,
             donor_meta = donor_meta, design_matrix = X, cohort = opt$cohort),
        file.path(opt$out_dir, sprintf("expression_INT_residuals_%s.rds", opt$cohort)))

# ============================================================================
# GENOTYPE DOSAGE MATRIX (from Schloissnig carrier lists in Supp Table 3 or 6)
# ============================================================================
st3 <- read_tsv(catalog_path, show_col_types = FALSE)
cat(sprintf("[%s] Loaded %s catalog: %d variants from %s\n",
            opt$cohort, opt$variant_class, nrow(st3), catalog_path))

classify_carriers <- function(row, cohort_set) {
  homs <- if (is.na(row$carriers_hom)) character(0) else strsplit(row$carriers_hom, ",")[[1]]
  hets <- if (is.na(row$carriers_het)) character(0) else strsplit(row$carriers_het, ",")[[1]]
  n_hom <- sum(homs %in% cohort_set); n_het <- sum(hets %in% cohort_set)
  af    <- (2*n_hom + n_het) / (2*length(cohort_set))
  data.frame(n_hom, n_het, n_ref = length(cohort_set) - n_hom - n_het,
             maf = pmin(af, 1 - af))
}

cohort_donors <- colnames(int_mat)
nz <- st3 %>% filter(trimws(NA12878_GT) == "0/0")
nz <- bind_cols(nz, bind_rows(lapply(seq_len(nrow(nz)),
                                     function(i) classify_carriers(nz[i,], cohort_donors))))
testable <- nz %>% filter(n_hom >= opt$min_hom, n_ref >= opt$min_ref, maf > opt$maf_floor)
cat(sprintf("\n[%s] Testable variants: %d  (NA12878=0/0 + MAF>%.2f + >=%d hom + >=%d ref)\n",
            opt$cohort, nrow(testable), opt$maf_floor, opt$min_hom, opt$min_ref))

n_var <- nrow(testable); n_donor <- length(cohort_donors)
dose_mat <- matrix(0L, n_var, n_donor, dimnames = list(testable$variant_id, cohort_donors))
for (i in seq_len(n_var)) {
  homs <- if (is.na(testable$carriers_hom[i])) character(0) else strsplit(testable$carriers_hom[i], ",")[[1]]
  hets <- if (is.na(testable$carriers_het[i])) character(0) else strsplit(testable$carriers_het[i], ",")[[1]]
  hom_idx <- match(intersect(homs, cohort_donors), cohort_donors)
  het_idx <- match(intersect(hets, cohort_donors), cohort_donors)
  if (length(hom_idx) > 0) dose_mat[i, hom_idx] <- 2L
  if (length(het_idx) > 0) dose_mat[i, het_idx] <- 1L
}
saveRDS(list(dose_mat = dose_mat, testable = testable, cohort = opt$cohort),
        file.path(opt$out_dir, sprintf("genotype_dosage_%s.rds", opt$cohort)))

# ============================================================================
# MATRIXEQTL
# ============================================================================
snpspos  <- data.frame(snp = testable$variant_id, chr = testable$chrom, pos = testable$pos)
gene_idx <- match(gene_unversioned, taylor$ensg_uv)
genepos  <- data.frame(geneid = gene_versioned,
                       chrom  = taylor$chrom[gene_idx],
                       s1     = taylor$tss[gene_idx],
                       s2     = taylor$tss[gene_idx])

snps <- SlicedData$new(); snps$CreateFromMatrix(dose_mat)
gene <- SlicedData$new(); gene$CreateFromMatrix(int_mat)
cvrt <- SlicedData$new()  # phenotype already covariate-residualized

out_cis <- tempfile(fileext = ".tsv")
me <- Matrix_eQTL_main(
  snps                  = snps, gene = gene, cvrt = cvrt,
  output_file_name      = NULL,           pvOutputThreshold     = 0,
  output_file_name.cis  = out_cis,        pvOutputThreshold.cis = 1,
  cisDist               = opt$cis_dist,   useModel              = modelLINEAR,
  errorCovariance       = numeric(),
  snpspos = snpspos, genepos = genepos,
  verbose = FALSE, pvalue.hist = FALSE)
cat(sprintf("[%s] MatrixEQTL: %d cis tests\n", opt$cohort, me$cis$ntests))

# ---- Annotate output + BH-FDR ----
cis_df <- read_tsv(out_cis, show_col_types = FALSE) %>%
  rename(variant_id = SNP, gene_id = gene, t_stat = `t-stat`,
         beta = `beta`, p = `p-value`, fdr_matrixeqtl = FDR)
cis_df$gene_id_uv <- sub("\\.[0-9]+$", "", cis_df$gene_id)
cis_df$q          <- p.adjust(cis_df$p, method = "BH")

gene_meta <- data.frame(gene_id = gene_versioned, gene_id_uv = gene_unversioned,
                        gene_symbol = taylor$geneSymbol[gene_idx],
                        gene_chrom  = taylor$chrom[gene_idx],
                        gene_TSS    = taylor$tss[gene_idx],
                        gene_strand = taylor$strand[gene_idx])
var_pos_lookup <- setNames(testable$pos, testable$variant_id)
cis_df <- cis_df %>%
  left_join(gene_meta, by = c("gene_id","gene_id_uv")) %>%
  mutate(variant_pos = var_pos_lookup[variant_id],
         distance    = variant_pos - gene_TSS)

# Cohort-specific carrier columns
cohort_suffix <- if (opt$cohort == "MAGE") "MAGE260" else "GEU121"
testable_short <- testable %>% select(variant_id, family, SVLEN,
                                      n_hom_x = n_hom, n_het_x = n_het, n_ref_x = n_ref, maf_x = maf)
colnames(testable_short) <- gsub("_x$", paste0("_", cohort_suffix), colnames(testable_short))
cis_df <- cis_df %>% left_join(testable_short, by = "variant_id")

out_tsv <- file.path(opt$out_dir, "all_polymorphicTE_eqtls.tsv")
write_tsv(cis_df, out_tsv)
cat(sprintf("\n[%s] Saved %s (%d rows)\n", opt$cohort, out_tsv, nrow(cis_df)))
cat(sprintf("[%s] Significance: q<0.05=%d  q<0.10=%d  q<0.20=%d\n",
            opt$cohort, sum(cis_df$q < 0.05), sum(cis_df$q < 0.10), sum(cis_df$q < 0.20)))
