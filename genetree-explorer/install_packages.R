# Install packages needed to run GeneTreeExplorer.
# Run once in R:  source("install_packages.R")
#
# CRAN packages come from the cloud CRAN mirror. Bioconductor packages
# (Biostrings, DECIPHER, ggtree) are installed with BiocManager.
# The README pins a Bioconductor release snapshot so lab machines stay in sync.

options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Pin a Bioconductor release that matches this R. 3.18 is the last release
# for R 4.3; newer R versions will take 3.19+. README documents the pin.
r_maj <- as.numeric(R.version$major)
r_min <- as.numeric(strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1])
bioc_version <- if (r_maj > 4 || (r_maj == 4 && r_min >= 5)) {
  "3.21"
} else if (r_maj == 4 && r_min >= 4) {
  "3.19"
} else {
  "3.18"
}
tryCatch(
  BiocManager::install(version = bioc_version, ask = FALSE, update = FALSE),
  error = function(e) {
    message(
      "Could not pin Bioconductor ", bioc_version, " (", conditionMessage(e),
      "). Using the default Bioconductor version for this R."
    )
  }
)

pkgs <- c(
  # User interface
  "shiny", "bslib", "DT", "shinycssloaders",
  # NCBI HTTP + JSON
  "httr2", "jsonlite", "memoise",
  # Data wrangling
  "dplyr", "purrr", "tibble",
  # Sequences, alignment, trees
  "Biostrings", "DECIPHER", "ape", "phangorn", "ggtree", "ggplot2",
  # Tests
  "testthat"
)

BiocManager::install(pkgs, ask = FALSE, update = FALSE, Ncpus = max(1L, parallel::detectCores(logical = FALSE) - 1L))

# Bioconductor 3.18 treeio can fail against current CRAN tidytree.
# Retry ggtree/treeio from GitHub in that case (Windows binary Bioc installs
# of a matching release usually do not need this).
if (!requireNamespace("ggtree", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  try(remotes::install_github("YuLab-SMU/treeio", upgrade = "never"), silent = TRUE)
  try(remotes::install_github("YuLab-SMU/ggtree", upgrade = "never"), silent = TRUE)
}

missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("These packages did not install: ", paste(missing, collapse = ", "))
}

message("All GeneTreeExplorer packages are installed.")
message("Bioconductor version: ", as.character(BiocManager::version()))
