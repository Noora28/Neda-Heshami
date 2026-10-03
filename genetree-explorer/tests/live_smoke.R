#!/usr/bin/env Rscript
# Optional live NCBI smoke test. Never fail the build if NCBI is down.
# Default path from the spec: pax6, Xenopus laevis, default species,
# protein, NJ, no bootstrap.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg)) {
  tests_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  tests_dir <- normalizePath("tests")
}
app_dir <- dirname(tests_dir)
setwd(app_dir)
source(file.path(app_dir, "R", "ncbi_api.R"), local = FALSE)
source(file.path(app_dir, "R", "phylo.R"), local = FALSE)

smoke_one <- function(symbol, taxon = "8355") {
  cat("SMOKE", symbol, "taxon", taxon, "\n")
  genes <- resolve_reference_genes(symbol, taxon)
  cat("  genes:", nrow(genes), paste(genes$symbol, collapse = ", "), "\n")
  ortho <- collect_ortholog_table(
    genes,
    comparison_taxids = unname(DEFAULT_SPECIES),
    use_all_orthologs = FALSE,
    max_n = MAX_ORTHOLOGS_DEFAULT
  )
  cat("  orthologs:", nrow(ortho), "\n")
  iso <- fetch_product_reports(ortho$gene_id)
  chosen <- choose_representative_isoform(iso, seq_type = "protein")
  chosen <- dplyr::left_join(
    ortho,
    chosen[, c("gene_id", "protein_acc", "protein_length", "transcript_acc")],
    by = "gene_id"
  )
  fasta <- efetch_fasta(chosen$protein_acc, db = "protein")
  seqs <- sequences_from_fasta_tbl(fasta, chosen, "protein")
  seqs <- seqs[!is.na(seqs$sequence) & nzchar(seqs$sequence), ]
  cat("  sequences:", nrow(seqs), "\n")
  if (nrow(seqs) < 3) {
    stop("fewer than 3 sequences")
  }
  aln <- align_sequences(seqs, seq_type = "protein")
  built <- build_tree(aln, seq_type = "protein", method = "NJ", bootstrap = 0)
  if (is.null(built$tree) || length(built$tree$tip.label) < 3) {
    stop("tree has fewer than 3 tips")
  }
  cat("  tree tips:", length(built$tree$tip.label), "\n")
  invisible(list(gene = genes, n_tips = length(built$tree$tip.label)))
}

run_safe <- function(symbol) {
  tryCatch(
    {
      smoke_one(symbol)
      list(symbol = symbol, ok = TRUE, detail = "ok")
    },
    error = function(e) {
      cat("  FAILED:", conditionMessage(e), "\n")
      list(symbol = symbol, ok = FALSE, detail = conditionMessage(e))
    }
  )
}

results <- list(run_safe("pax6"))
if (isTRUE(results[[1]]$ok)) {
  results <- c(results, list(run_safe("sox2"), run_safe("shh")))
}

out <- file.path(app_dir, "tests", "live_smoke_result.txt")
lines <- vapply(results, function(r) {
  sprintf("%s\t%s\t%s", r$symbol, if (isTRUE(r$ok)) "PASS" else "FAIL", r$detail)
}, character(1))
stamp <- format(Sys.time(), tz = "UTC", usetz = TRUE)
writeLines(c(paste("recorded", stamp), lines), out)
cat(paste(c(paste("recorded", stamp), lines), collapse = "\n"), "\n")
# Always exit 0: live NCBI must not block.
quit(status = 0)
