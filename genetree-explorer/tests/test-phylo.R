# Alignment / NJ tests using tiny in-memory sequences (no NCBI).

testthat::test_that("n < 3 skips the tree and explains why", {
  reason <- skip_tree_reason(2)
  testthat::expect_match(reason, "at least 3")
  seqs <- Biostrings::AAStringSet(c(a = "MKTAYIAK", b = "MKTAYIAQ"))
  aln <- DECIPHER::AlignSeqs(seqs, verbose = FALSE, processors = 1)
  built <- build_tree(aln, seq_type = "protein", method = "NJ")
  testthat::expect_null(built$tree)
  testthat::expect_match(built$message, "at least 3")
})

testthat::test_that("X. laevis L and S remain separate unique tips", {
  genes <- tibble::tibble(
    gene_id = c("379100", "100337586", "5080"),
    symbol = c("pax6.L", "pax6.S", "PAX6"),
    taxname = c("Xenopus laevis", "Xenopus laevis", "Homo sapiens")
  )
  labs <- make_tip_labels(genes)
  testthat::expect_equal(length(labs), 3)
  testthat::expect_equal(length(unique(labs)), 3)
  testthat::expect_true(any(grepl("pax6\\.L", labs)))
  testthat::expect_true(any(grepl("pax6\\.S", labs)))
})

testthat::test_that("NJ protein tree has at least 3 tips", {
  seqs <- Biostrings::AAStringSet(c(
    pax6.L_Xenopus_laevis = "MKTAYIAKQRQISFVKSHFS",
    pax6.S_Xenopus_laevis = "MKTAYIAKQRQISFVKSHFA",
    PAX6_Homo_sapiens = "MKTAYVAKQRQISFVKSHFS",
    Pax6_Mus_musculus = "MKTAYIAKQRQIAFVKSHFS"
  ))
  aln <- DECIPHER::AlignSeqs(seqs, verbose = FALSE, processors = 1)
  built <- build_tree(aln, seq_type = "protein", method = "NJ", bootstrap = 0)
  testthat::expect_null(built$message)
  testthat::expect_s3_class(built$tree, "phylo")
  testthat::expect_gte(length(built$tree$tip.label), 3)
  ident <- mean_pairwise_identity(aln)
  testthat::expect_true(ident > 0.5 && ident <= 1)
})

testthat::test_that("CDS without a coding frame is refused", {
  cds <- Biostrings::DNAStringSet(c(
    a = "ATGAAATAG",
    b = "ATGAAAA"
  ))
  testthat::expect_error(translate_cds_or_stop(cds), "coding frame")
})

testthat::test_that("ML is refused when n > 12", {
  testthat::expect_match(ml_refused_reason(13), "refused")
  testthat::expect_null(ml_refused_reason(5))
})

testthat::test_that("UPGMA returns a tree for 4 sequences", {
  seqs <- Biostrings::AAStringSet(c(
    t1 = "MKTAYIAKQRQISFVKSHFS",
    t2 = "MKTAYIAKQRQISFVKSHFA",
    t3 = "MKTAYVAKQRQISFVKSHFS",
    t4 = "MKTAYIAKQRQIAFVKSHFS"
  ))
  aln <- DECIPHER::AlignSeqs(seqs, verbose = FALSE, processors = 1)
  built <- build_tree(aln, seq_type = "protein", method = "UPGMA")
  testthat::expect_s3_class(built$tree, "phylo")
})
