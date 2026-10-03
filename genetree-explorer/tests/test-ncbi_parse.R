# Fixture tests for NCBI JSON parsers, isoform choice, homeologs, and caps.
# Run from genetree-explorer/:  Rscript tests/run_tests.R

.fixture <- function(name) {
  file.path(.test_dir, "fixtures", name)
}

read_json_fixture <- function(name) {
  jsonlite::fromJSON(.fixture(name), simplifyVector = FALSE)
}

testthat::test_that("Xenopus pax6 symbol lookup returns L and S homeologs", {
  json <- read_json_fixture("gene_symbol_pax6_xenopus.json")
  tbl <- parse_gene_reports(json)
  testthat::expect_true(nrow(tbl) >= 2)
  testthat::expect_true(all(c("pax6.L", "pax6.S") %in% tbl$symbol))
  testthat::expect_true(all(tbl$tax_id == "8355"))
  testthat::expect_true("379100" %in% tbl$gene_id)
  testthat::expect_true("100337586" %in% tbl$gene_id)
})

testthat::test_that("gene-by-id fixture parses human PAX6", {
  json <- read_json_fixture("gene_id_5080.json")
  tbl <- parse_gene_reports(json)
  testthat::expect_equal(nrow(tbl), 1)
  testthat::expect_equal(tbl$gene_id[[1]], "5080")
  testthat::expect_equal(tbl$symbol[[1]], "PAX6")
  testthat::expect_equal(tbl$taxname[[1]], "Homo sapiens")
})

testthat::test_that("missing gene fields become NA instead of crashing", {
  tbl <- parse_gene_reports(list(reports = list(list(gene = list()))))
  testthat::expect_equal(nrow(tbl), 1)
  testthat::expect_true(is.na(tbl$gene_id[[1]]) || !nzchar(tbl$gene_id[[1]]))
})

testthat::test_that("taxon_suggest parses sci_name_and_ids", {
  json <- read_json_fixture("taxon_suggest_xenopus.json")
  tbl <- parse_taxon_suggest(json)
  testthat::expect_true(nrow(tbl) >= 1)
  testthat::expect_true("Xenopus laevis" %in% tbl$sci_name)
  testthat::expect_true("8355" %in% tbl$tax_id)
})

testthat::test_that("esearch JSON yields gene id 5080", {
  json <- read_json_fixture("esearch_pax6_human.json")
  ids <- parse_esearch_ids(json)
  testthat::expect_equal(ids, "5080")
})

testthat::test_that("isoform choice prefers NP_ over a longer XP_", {
  json <- read_json_fixture("product_report_379100.json")
  iso <- parse_product_isoforms(json)
  testthat::expect_true(any(grepl("^NP_", iso$protein_acc)))
  testthat::expect_true(any(grepl("^XP_", iso$protein_acc)))
  chosen <- choose_representative_isoform(iso, seq_type = "protein")
  testthat::expect_equal(nrow(chosen), 1)
  testthat::expect_match(chosen$protein_acc[[1]], "^NP_")
  # The model XP_ isoform in the fixture is longer (467 vs 453).
  xp_len <- max(iso$protein_length[grepl("^XP_", iso$protein_acc)], na.rm = TRUE)
  np_len <- chosen$protein_length[[1]]
  testthat::expect_true(xp_len > np_len)
})

testthat::test_that("CDS isoform choice prefers NM_ over XM_", {
  json <- read_json_fixture("product_report_379100.json")
  iso <- parse_product_isoforms(json)
  chosen <- choose_representative_isoform(iso, seq_type = "cds")
  testthat::expect_match(chosen$transcript_acc[[1]], "^NM_")
})

testthat::test_that("human PAX6 representative is a curated NP_ accession", {
  json <- read_json_fixture("product_report_5080.json")
  iso <- parse_product_isoforms(json)
  chosen <- choose_representative_isoform(iso, seq_type = "protein")
  testthat::expect_match(chosen$protein_acc[[1]], "^NP_")
  testthat::expect_equal(chosen$protein_length[[1]], max(iso$protein_length, na.rm = TRUE))
})

testthat::test_that("ortholog cap keeps reference and X. laevis L/S first", {
  json <- read_json_fixture("orthologs_5080_page.json")
  tbl <- parse_gene_reports(json)
  xl <- tibble::tibble(
    gene_id = c("379100", "100337586"),
    symbol = c("pax6.L", "pax6.S"),
    description = "paired box 6",
    tax_id = "8355",
    taxname = "Xenopus laevis",
    common_name = "African clawed frog",
    type = "PROTEIN_CODING",
    chromosomes = c("4L", "4S"),
    summary = NA_character_,
    ensembl_gene_ids = NA_character_,
    swiss_prot = NA_character_
  )
  merged <- dplyr::bind_rows(xl, tbl)
  capped <- cap_orthologs(
    merged,
    reference_gene_ids = c("379100", "100337586"),
    max_n = 6
  )
  testthat::expect_equal(nrow(capped), 6)
  testthat::expect_equal(capped$symbol[1:2], c("pax6.L", "pax6.S"))
})

testthat::test_that("parse_fasta_text reads headers and concatenates lines", {
  txt <- paste(
    ">NP_000271.1 paired box [Homo sapiens]",
    "MKTAA",
    "LLAAA",
    ">XP_FAKE.1 model",
    "MKTAV",
    sep = "\n"
  )
  fa <- parse_fasta_text(txt)
  testthat::expect_equal(nrow(fa), 2)
  testthat::expect_equal(fa$accession[[1]], "NP_000271.1")
  testthat::expect_equal(fa$sequence[[1]], "MKTAALLAAA")
})
