#!/usr/bin/env Rscript
# Install the R packages used by the eQTL pipeline + R figure/supp-table scripts.
# Idempotent: skips already-installed packages.
#
# Example usage:
# Rscript scripts/install_R_deps.R

required <- c(
  # core eQTL pipeline
  "MatrixEQTL", "matrixStats", "R.utils", "optparse", "data.table",
  "dplyr", "readr",
  # plotting (figures + supp figures)
  "ggplot2", "ggrepel", "ggtext", "patchwork", "cowplot", "tidyr", "tibble"
)
to_install <- setdiff(required, rownames(installed.packages()))
if (length(to_install) == 0) {
  cat("All required packages already installed:\n")
  for (p in required) cat("  ", p, " ", as.character(packageVersion(p)), "\n", sep = "")
  quit(status = 0)
}
cat("Installing:", paste(to_install, collapse = ", "), "\n")
install.packages(to_install, repos = "https://cloud.r-project.org")
for (p in required) {
  if (!requireNamespace(p, quietly = TRUE)) stop("FAILED to install ", p)
  cat("  ", p, " ", as.character(packageVersion(p)), "\n", sep = "")
}
cat("Done.\n")
