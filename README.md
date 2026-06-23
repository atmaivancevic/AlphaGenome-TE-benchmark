# AlphaGenome-TE-benchmark

Scripts and workflows used in "Benchmarking AlphaGenome on repeat-derived regulatory variants across the human genome" (2026) Atma Ivancevic, Daniel Larremore, Ryan Layer, Edward B. Chuong

---

### Programs Used

- **AlphaGenome** v0.6.1 https://github.com/google-deepmind/alphagenome
- **Python** v3.11 https://www.python.org/
- **R** v4.6.0 https://www.r-project.org/
- **MatrixEQTL** v2.3 https://github.com/andreyshabalin/MatrixEQTL
- **bcftools** v1.23.1 http://www.htslib.org/
- **bedtools** v2.31.1 https://github.com/arq5x/bedtools2
- **Samtools** v1.16.1 http://www.htslib.org/
- **GIGGLE** v0.6.3 https://github.com/ryanlayer/giggle

Python packages are pinned in [requirements.txt](requirements.txt); R packages install via [install_R_deps.R](scripts/install_R_deps.R). Each script's header has an `Example usage` block with the exact command used.

---

### Setup

```bash
python -m venv venv && source venv/bin/activate
pip install -r requirements.txt
Rscript scripts/install_R_deps.R

# AlphaGenome API key (required for all prediction/scoring scripts)
cp scripts/my_api_key.txt.example scripts/my_api_key.txt   # paste your key
```

Run all scripts from the repository root.

---

### Fig 1 — Genome-wide H3K27ac prediction across ENCODE biosamples

1. **Build the peak manifest, download, and merge**
   1) [resolve_encode_h3k27ac_peaks.sh](scripts/fig1_baseline_chromatin/resolve_encode_h3k27ac_peaks.sh) — map each AG-training bigWig to its ENCODE narrowPeak
   2) [fetch_encode_h3k27ac_peaks.sh](scripts/fig1_baseline_chromatin/fetch_encode_h3k27ac_peaks.sh) — download + md5-verify peak files
   3) [merge_encode_h3k27ac_peaks.sh](scripts/fig1_baseline_chromatin/merge_encode_h3k27ac_peaks.sh) — merge multi-experiment biosamples (calls [group_peaks_by_biosample.py](scripts/fig1_baseline_chromatin/group_peaks_by_biosample.py))
   4) [fetch_encode_h3k27ac_bigwigs.sh](scripts/fig1_baseline_chromatin/fetch_encode_h3k27ac_bigwigs.sh) — fetch per-biosample fold-change bigWigs

2. **AG prediction + correlation (Panel A)**
   1) [predict_peak_signal_batched.py](scripts/fig1_baseline_chromatin/predict_peak_signal_batched.py) — predict AG H3K27ac at each experimental peak
   2) [plot_peak_correlation.py](scripts/fig1_baseline_chromatin/plot_peak_correlation.py) — AG-vs-experimental scatter (Pearson r + Spearman ρ)

3. **GIGGLE TE-family enrichment (Panel D + Supp Fig S2)**
   1) [giggle_enrichment.sh](scripts/fig1_baseline_chromatin/giggle_enrichment.sh) — GIGGLE vs. the hg38 RepeatMasker index (run on HPC)
   2) [plot_giggle_bubbles.py](scripts/fig1_baseline_chromatin/plot_giggle_bubbles.py) — bubble plots

4. **Panel B/C outliers + browser shots**
   1) [find_panel_outliers.py](scripts/fig1_baseline_chromatin/find_panel_outliers.py) — rank the most over-/under-predicted peaks
   2) [annotate_panel_outliers.py](scripts/fig1_baseline_chromatin/annotate_panel_outliers.py) — nearest gene, TSS distance, promoter/proximal/distal class
   3) [plot_panel_BC_browser.py](scripts/fig1_baseline_chromatin/plot_panel_BC_browser.py) — 3-track browser shots (NFKBIA over-prediction, PCDH7 under-prediction)

---

### Fig 2 — Polymorphic TE insertion example (LTR5_Hs at HLA-DQA2)

1. **Build the variant and run AlphaGenome**
   1) [filter_solo_TE_INS_vcf.py](scripts/fig2_polymorphic_TE_example/filter_solo_TE_INS_vcf.py) — filter the SVAN-annotated BCF to solo MEI insertions
   2) [fetch_variant_alt_seq.py](scripts/fig2_polymorphic_TE_example/fetch_variant_alt_seq.py) — pull the anchor + ALT sequence for the exemplar variant
   3) [plot_polymorphic_TE_insertion.py](scripts/fig2_polymorphic_TE_example/plot_polymorphic_TE_insertion.py) — AG-predicted RNA/ChIP/ATAC tracks (ref vs. alt) at the insertion

---

### Fig 3 — Polymorphic TE cis-eQTLs (MAGE / GEUVADIS)

1. **Run the cis-eQTL pipeline**
   1) [eqtl_matrixeqtl_pipeline.R](scripts/fig3_polymorphic_TE_eQTLs/eqtl_matrixeqtl_pipeline.R) — MatrixEQTL cis-eQTL (MAGE-260 / GEUVADIS-121 × INS/DEL)
   2) [cross_cohort_compare.R](scripts/fig3_polymorphic_TE_eQTLs/cross_cohort_compare.R) — MAGE-vs-GEUVADIS side-by-side table

2. **Score variants with AlphaGenome and merge**
   1) [score_variant_lfc.py](scripts/fig4_5_LTR10_CRISPR_comparison/score_variant_lfc.py) — AG RNA log2FC scoring
   2) [score_variant_chromatin.py](scripts/fig3_polymorphic_TE_eQTLs/score_variant_chromatin.py) — AG chromatin (H3K27ac/H3K4me1/ATAC) scoring
   3) [merge_AG_scores_into_supp_tables.py](scripts/fig3_polymorphic_TE_eQTLs/merge_AG_scores_into_supp_tables.py) — merge AG scores into Supp Tables 4/5/7/8

3. **Plots**
   1) [plot_fig3_eqtl_boxplots.R](scripts/fig3_polymorphic_TE_eQTLs/plot_fig3_eqtl_boxplots.R) — per-candidate eQTL boxplots
   2) [plot_fig3_scatter_AG_vs_MAGE.R](scripts/fig3_polymorphic_TE_eQTLs/plot_fig3_scatter_AG_vs_MAGE.R) — AG-predicted vs. MAGE-observed effect-size scatter
   3) [render_fig3_AG_panels.sh](scripts/fig3_polymorphic_TE_eQTLs/render_fig3_AG_panels.sh) — batch-render AG genome-browser panels for the candidates
   4) [compose_fig3.R](scripts/fig3_polymorphic_TE_eQTLs/compose_fig3.R) — assemble the final multi-panel Fig 3
   5) [plot_supp_fig_AG_per_track.R](scripts/fig3_polymorphic_TE_eQTLs/plot_supp_fig_AG_per_track.R) — per-RNA-track AG supplementary figure

---

### Fig 4 & 5 — LTR10 CRISPRi enhancer benchmarking (HCT116)

1. **Score the CRISPRi-validated LTR10 enhancers**
   1) [score_variant_lfc.py](scripts/fig4_5_LTR10_CRISPR_comparison/score_variant_lfc.py) — AG RNA log2FC scoring
   2) [score_variant_chromatin.py](scripts/fig3_polymorphic_TE_eQTLs/score_variant_chromatin.py) — AG chromatin scoring

2. **Plots**
   1) [plot_LTR10AF_chromatin_waterfall.py](scripts/fig4_5_LTR10_CRISPR_comparison/plot_LTR10AF_chromatin_waterfall.py) — 650-element chromatin waterfall (Fig 5A)
   2) [plot_LTR10AF_experimental_vs_AG_scatter.py](scripts/fig4_5_LTR10_CRISPR_comparison/plot_LTR10AF_experimental_vs_AG_scatter.py) — experimental vs. AG ChIP scatter (Fig 5B)
   3) [plot_fig5_scatter_AG_vs_CRISPR.R](scripts/fig4_5_LTR10_CRISPR_comparison/plot_fig5_scatter_AG_vs_CRISPR.R) — per-enhancer AG-vs-CRISPRi effect-size scatters
   4) [plot_LTR10_dumbbell.R](scripts/fig4_5_LTR10_CRISPR_comparison/plot_LTR10_dumbbell.R) — CRISPRi-vs-AG dumbbell per enhancer
   5) [plot_AG_browser_shot.py](scripts/fig4_5_LTR10_CRISPR_comparison/plot_AG_browser_shot.py) — AG genome-browser track shot
   6) [recolor_AG_screenshot_tracks.py](scripts/fig4_5_LTR10_CRISPR_comparison/recolor_AG_screenshot_tracks.py) — recolour browser-shot tracks to the locked assay palette
   7) [plot_variant_coord_scatter.py](scripts/fig4_5_LTR10_CRISPR_comparison/plot_variant_coord_scatter.py) — AG effect vs. gene-TSS coordinate scatter

---

### Fig 6 — In-silico AP1 motif perturbation at LTR10.ATG12

1. **Generate perturbation alleles** (deterministic, seed-locked)
   1) [generate_LTR10_ATG12_AP1_perturbations.py](scripts/fig6_AP1_perturbation/generate_LTR10_ATG12_AP1_perturbations.py) — scramble + TF-substitution alleles
   2) [generate_LTR10_ATG12_AP1_insertions.py](scripts/fig6_AP1_perturbation/generate_LTR10_ATG12_AP1_insertions.py) — motif-addition alleles (past WT density)

2. **Run AlphaGenome**
   1) [predict_variant_tracks.py](scripts/fig6_AP1_perturbation/predict_variant_tracks.py) — per-allele AG signal (TF / chromatin / RNA)

3. **Plots**
   1) [plot_fig6_panelA1_TF_titration.py](scripts/fig6_AP1_perturbation/plot_fig6_panelA1_TF_titration.py) — TF binding titration (Panel A1)
   2) [plot_fig6_panelA2_chromatin_RNA_titration.py](scripts/fig6_AP1_perturbation/plot_fig6_panelA2_chromatin_RNA_titration.py) — chromatin + RNA titration (Panel A2)
   3) [plot_fig6_panelB_TF_specificity.py](scripts/fig6_AP1_perturbation/plot_fig6_panelB_TF_specificity.py) — TF-substitution specificity bars (Panel B)

---

### Supplementary tables

- [build_supp_table1_encode_h3k27ac_tracks.py](scripts/generate_supp_tables/build_supp_table1_encode_h3k27ac_tracks.py) — Supp 1: ENCODE H3K27ac track manifest
- [build_supp_table2_individuals.py](scripts/generate_supp_tables/build_supp_table2_individuals.py) — Supp 2: per-donor cohort flags (908 individuals)
- [build_supp_table3_polymorphic_TE_insertions.py](scripts/generate_supp_tables/build_supp_table3_polymorphic_TE_insertions.py) — Supp 3: polymorphic TE insertion catalog
- [build_supp_table6_polymorphic_TE_deletions.py](scripts/generate_supp_tables/build_supp_table6_polymorphic_TE_deletions.py) — Supp 6: polymorphic TE deletion catalog
- [build_supp_tables_4_5_7_8_eqtl.R](scripts/generate_supp_tables/build_supp_tables_4_5_7_8_eqtl.R) — Supp 4/5 (INS) + 7/8 (DEL): eQTL summary tables
- [build_supp_table10_LTR10_CRISPR_AG_predictions.py](scripts/generate_supp_tables/build_supp_table10_LTR10_CRISPR_AG_predictions.py) — Supp 10: LTR10 CRISPRi vs. AG predictions
- [build_supp_table11_LTR10AF_experimental_vs_AG.py](scripts/generate_supp_tables/build_supp_table11_LTR10AF_experimental_vs_AG.py) — Supp 11: 650 LTR10A/F experimental vs. AG

(Supp Table 9, the eQTL PC × covariate diagnostic, is written by `eqtl_matrixeqtl_pipeline.R`.)

---

### Notes

- All software was run with default settings unless otherwise indicated; each script's header gives the exact command used to produce the paper results.
- Run on macOS; the GIGGLE enrichment and bulk ENCODE/ENA downloads were run on a Linux HPC.
- Several plotting and table scripts read the supplementary-table TSVs (`supptables/`) or curated input tables (`data/`) as inputs. Run the relevant `generate_supp_tables/` builder first, or place the published supplementary TSVs in `supptables/`.
- eQTL caveats: the analysis is restricted to NA12878 = 0/0 variants, to avoid AG's GM12878 reference-allele bias. MAGE and GEUVADIS share some 1000G donors, so GEUVADIS is technical (not independent) replication. See the manuscript Methods.

---

### Data

- **ENCODE H3K27ac** peaks + bigWigs — downloaded via the URLs in Supp Table 1 (built by `build_supp_table1_encode_h3k27ac_tracks.py`).
- **Schloissnig 1KG-ONT-Vienna polymorphic TEs** (IGSR `1KG_ONT_VIENNA` release v1.1):
  - `svim.asm.hg38.noGt.SVAN_1.3.bcf` — md5 `fe816142d21ef48960fd9f6ecb82ec6e`
  - `svim.asm.hg38.bcf` (908 genotypes) — md5 `ac4751b69374dee1cec8414d771947ee`
  - ShapeIt5-phased callset (908-sample list) — md5 `56361b82acb52f96ef5f176b4c6c3fad`
  - Base URL: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/release/v1.1/
- **Schloissnig sample metadata** (`sample.tsv`) — Zenodo https://zenodo.org/records/14535469 (`sv-analysis-v1.1.zip`, md5 `a36f30379d34b29a302517a39d568e59`)
- **MAGE** expression (Taylor et al. 2024) — https://github.com/mccoy-lab/MAGE (Zenodo [10.5281/zenodo.10535719](https://zenodo.org/records/10535719))
- **GEUVADIS** RPKM — https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEUV/E-GEUV-1/analysis_results/GD660.GeneQuantRPKM.txt.gz
- **GENCODE v46** annotation (Feather) — auto-downloaded by the AG plotting scripts on first run.

---

### License

MIT License — see [LICENSE](LICENSE).
