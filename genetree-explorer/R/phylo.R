# Alignment and tree helpers (no Shiny dependency).
# v1 path that should always work when NCBI returned enough sequences:
#   protein sequences -> DECIPHER::AlignSeqs -> neighbour-joining (NJ)
#   -> midpoint root (unless a unique outgroup tip is named).
#
# Extra methods (UPGMA, ML) are available but guarded so they cannot hang
# the app: ML is refused for n > 12, and ML never runs bootstrap.

if (!exists("MAX_ORTHOLOGS_HARD", inherits = TRUE)) {
  source(file.path(dirname(normalizePath(sys.frame(1)$ofile)), "config.R"), local = FALSE)
}

sanitize_tip_label <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9._-]", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x[nchar(x) == 0 | is.na(x)] <- "tip"
  x
}

#' Build unique tip labels. Keep X. laevis .L and .S distinct.
make_tip_labels <- function(genes) {
  raw <- paste(genes$symbol, gsub(" ", "_", genes$taxname), sep = "_")
  raw <- sanitize_tip_label(raw)
  # If two genes share symbol+species (homeologs already differ by .L/.S),
  # still guarantee uniqueness with GeneID.
  dups <- duplicated(raw) | duplicated(raw, fromLast = TRUE)
  raw[dups] <- sanitize_tip_label(paste(raw[dups], genes$gene_id[dups], sep = "_"))
  make.unique(raw, sep = "_")
}

skip_tree_reason <- function(n) {
  n <- as.integer(n)
  if (is.na(n) || n < 3) {
    paste0(
      "Need at least 3 sequences to build a tree (have ",
      ifelse(is.na(n), 0, n),
      "). Fetch more species or choose a gene with more NCBI orthologs."
    )
  } else {
    NULL
  }
}

ml_refused_reason <- function(n) {
  if (as.integer(n) > ML_MAX_SEQUENCES) {
    paste0(
      "Maximum-likelihood is refused when more than ", ML_MAX_SEQUENCES,
      " sequences are selected (have ", n, "). Use NJ, or lower the ortholog cap."
    )
  } else {
    NULL
  }
}

#' Mean pairwise identity over columns where both residues are not gaps.
mean_pairwise_identity <- function(aln) {
  if (is.null(aln) || length(aln) < 2) return(NA_real_)
  mat <- as.matrix(aln)
  n <- nrow(mat)
  vals <- numeric(0)
  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      a <- mat[i, ]
      b <- mat[j, ]
      gap <- a %in% c("-", ".", "?") | b %in% c("-", ".", "?")
      use <- !gap
      if (!any(use)) next
      vals <- c(vals, mean(a[use] == b[use]))
    }
  }
  if (!length(vals)) NA_real_ else mean(vals)
}

alignment_to_fasta_text <- function(aln) {
  if (is.null(aln) || !length(aln)) return("")
  paste0(">", names(aln), "\n", as.character(aln), collapse = "\n")
}

unaligned_to_fasta_text <- function(seqs) {
  alignment_to_fasta_text(seqs)
}

#' Translate CDS; refuse if length is not a multiple of 3 (no coding frame).
translate_cds_or_stop <- function(cds_set) {
  widths <- Biostrings::width(cds_set)
  bad <- names(cds_set)[widths %% 3L != 0]
  if (length(bad)) {
    stop(
      "CDS translation failure: these sequences are not a complete coding frame (length not divisible by 3): ",
      paste(bad, collapse = ", "),
      ". Use protein sequences instead.",
      call. = FALSE
    )
  }
  Biostrings::translate(cds_set, if.fuzzy.codon = "solve")
}

#' Align proteins with DECIPHER, then emit codon alignment for the CDS.
protein_guided_cds_alignment <- function(cds_set) {
  aa <- translate_cds_or_stop(cds_set)
  aa_aln <- DECIPHER::AlignSeqs(aa, verbose = FALSE, processors = 1)
  # Back-translate: each aligned residue maps to one codon; gaps to ---.
  out <- character(length(cds_set))
  names(out) <- names(cds_set)
  for (nm in names(cds_set)) {
    prot_aln <- strsplit(as.character(aa_aln[[nm]]), "")[[1]]
    cds <- as.character(cds_set[[nm]])
    codon_i <- 1L
    pieces <- character(length(prot_aln))
    for (k in seq_along(prot_aln)) {
      if (prot_aln[[k]] == "-") {
        pieces[[k]] <- "---"
      } else {
        start <- (codon_i - 1L) * 3L + 1L
        pieces[[k]] <- substr(cds, start, start + 2L)
        codon_i <- codon_i + 1L
      }
    }
    out[[nm]] <- paste(pieces, collapse = "")
  }
  Biostrings::DNAStringSet(out)
}

align_sequences <- function(seq_tbl, seq_type = c("protein", "cds")) {
  seq_type <- match.arg(seq_type)
  seq_tbl <- seq_tbl[!is.na(seq_tbl$sequence) & nzchar(seq_tbl$sequence), , drop = FALSE]
  if (!nrow(seq_tbl)) {
    stop("No valid sequences were downloaded from NCBI.", call. = FALSE)
  }
  labels <- seq_tbl$tip_label
  seqs <- seq_tbl$sequence
  names(seqs) <- labels
  if (seq_type == "protein") {
    aset <- Biostrings::AAStringSet(seqs)
    aln <- DECIPHER::AlignSeqs(aset, verbose = FALSE, processors = 1)
  } else {
    aset <- Biostrings::DNAStringSet(seqs)
    aln <- protein_guided_cds_alignment(aset)
  }
  aln
}

to_phyDat <- function(aln, seq_type) {
  mat <- as.matrix(aln)
  rownames(mat) <- names(aln)
  if (identical(seq_type, "protein")) {
    phangorn::phyDat(mat, type = "AA")
  } else {
    # CDS alignment is nucleotide (codon-aware gaps of length 3). K80 needs DNA.
    phangorn::phyDat(mat, type = "DNA")
  }
}

dist_for <- function(pd, seq_type) {
  if (identical(seq_type, "protein")) {
    phangorn::dist.ml(pd, model = "JTT")
  } else {
    phangorn::dist.ml(pd, model = "K80")
  }
}

root_tree <- function(tree, outgroup_tip = NULL) {
  tips <- tree$tip.label
  if (!is.null(outgroup_tip) && nzchar(outgroup_tip) && outgroup_tip %in% tips) {
    if (sum(tips == outgroup_tip) != 1) {
      warning("Outgroup tip is not unique; using midpoint root.")
      return(phangorn::midpoint(tree))
    }
    ape::root(tree, outgroup = outgroup_tip, resolve.root = TRUE)
  } else {
    phangorn::midpoint(tree)
  }
}

#' Neighbour-joining (default v1 method).
build_nj_tree <- function(pd, seq_type, bootstrap = 0L, outgroup_tip = NULL) {
  dm <- dist_for(pd, seq_type)
  tree <- phangorn::NJ(dm)
  tree <- root_tree(tree, outgroup_tip)
  if (bootstrap > 0) {
    fun <- function(x) phangorn::NJ(dist_for(x, seq_type))
    bs <- phangorn::bootstrap.phyDat(pd, FUN = fun, bs = as.integer(bootstrap))
    tree <- phangorn::plotBS(tree, bs, type = "none", p = 0)
  }
  tree
}

#' UPGMA assumes a molecular clock (same rate on every branch).
build_upgma_tree <- function(pd, seq_type, bootstrap = 0L, outgroup_tip = NULL) {
  dm <- dist_for(pd, seq_type)
  tree <- phangorn::upgma(dm)
  tree <- root_tree(tree, outgroup_tip)
  if (bootstrap > 0) {
    fun <- function(x) phangorn::upgma(dist_for(x, seq_type))
    bs <- phangorn::bootstrap.phyDat(pd, FUN = fun, bs = as.integer(bootstrap))
    tree <- phangorn::plotBS(tree, bs, type = "none", p = 0)
  }
  tree
}

#' ML from an NJ start. No bootstrap (too slow; not offered in the UI).
build_ml_tree <- function(pd, seq_type, outgroup_tip = NULL) {
  n <- length(pd)
  reason <- ml_refused_reason(n)
  if (!is.null(reason)) stop(reason, call. = FALSE)
  dm <- dist_for(pd, seq_type)
  start <- phangorn::NJ(dm)
  model <- if (identical(seq_type, "protein")) "JTT" else "K80"
  fit <- phangorn::pml(start, data = pd, model = model)
  fit <- phangorn::optim.pml(
    fit,
    model = model,
    optNni = TRUE,
    rearrangement = "NNI",
    control = phangorn::pml.control(trace = 0)
  )
  tree <- root_tree(fit$tree, outgroup_tip)
  attr(tree, "pml_logLik") <- fit$logLik
  tree
}

build_tree <- function(aln, seq_type = "protein", method = c("NJ", "UPGMA", "ML"),
                       bootstrap = 0L, outgroup_tip = NULL) {
  method <- match.arg(method)
  n <- length(aln)
  reason <- skip_tree_reason(n)
  if (!is.null(reason)) {
    return(list(tree = NULL, message = reason, method = method))
  }
  if (identical(method, "ML") && bootstrap > 0) {
    bootstrap <- 0L
  }
  pd <- to_phyDat(aln, seq_type)
  tree <- switch(
    method,
    NJ = build_nj_tree(pd, seq_type, bootstrap = bootstrap, outgroup_tip = outgroup_tip),
    UPGMA = build_upgma_tree(pd, seq_type, bootstrap = bootstrap, outgroup_tip = outgroup_tip),
    ML = build_ml_tree(pd, seq_type, outgroup_tip = outgroup_tip)
  )
  list(tree = tree, message = NULL, method = method, n = n, bootstrap = as.integer(bootstrap))
}

plot_ggtree <- function(tree, layout = c("rectangular", "circular"),
                        show_bootstrap = TRUE) {
  layout <- match.arg(layout)
  if (is.null(tree)) {
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::ggtitle("No tree to display."))
  }
  p <- ggtree::ggtree(tree, layout = layout)
  p <- p + ggtree::geom_tiplab(size = 3.2, align = identical(layout, "rectangular"))
  # Bootstrap as node labels: plotBS stores them in tree$node.label.
  if (show_bootstrap && !identical(layout, "circular") &&
      !is.null(tree$node.label) && any(nzchar(as.character(tree$node.label)))) {
    p <- p + ggtree::geom_nodelab(size = 2.5, nudge_x = 0.01)
  }
  p + ggplot2::theme(plot.margin = ggplot2::margin(10, 80, 10, 10))
}

sequences_from_fasta_tbl <- function(fasta_tbl, genes, seq_type) {
  acc_col <- if (identical(seq_type, "protein")) "protein_acc" else "transcript_acc"
  # Match efetch accession (may include version) to chosen isoform.
  fasta_tbl$acc_core <- sub("\\..*$", "", fasta_tbl$accession)
  genes$acc_full <- genes[[acc_col]]
  genes$acc_core <- sub("\\..*$", "", genes$acc_full)
  genes$tip_label <- make_tip_labels(genes)
  seq_map <- fasta_tbl$sequence
  names(seq_map) <- fasta_tbl$accession
  seq_map2 <- fasta_tbl$sequence
  names(seq_map2) <- fasta_tbl$acc_core
  genes$sequence <- unname(seq_map[genes$acc_full])
  miss <- is.na(genes$sequence)
  genes$sequence[miss] <- unname(seq_map2[genes$acc_core[miss]])
  genes
}
