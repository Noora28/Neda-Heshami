# NCBI Datasets v2 + E-utilities helpers (no Shiny dependency).
#
# Confirmed against NCBI Datasets OpenAPI (github.com/ncbi/datasets
# datasets.openapi.yaml) and https://www.ncbi.nlm.nih.gov/books/NBK25499/
# (E-utilities, last update 4 March 2026):
#
#   Datasets base: https://api.ncbi.nlm.nih.gov/datasets/v2
#   GET /gene/symbol/{symbols}/taxon/{taxon}
#   GET /gene/id/{gene_ids}                    (still works; marked deprecated)
#   GET /gene/id/{gene_ids}/dataset_report     (current gene-report alias)
#   GET /gene/id/{gene_id}/orthologs
#   GET /gene/id/{gene_ids}/product_report     JSON: reports[].product
#   GET /taxonomy/taxon_suggest/{taxon_query}  JSON: sci_name_and_ids[]
#
#   E-utilities: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/
#   esearch.fcgi / esummary.fcgi  retmode=json
#   efetch.fcgi  db=protein|nuccore  rettype=fasta  retmode=text
#   Datasets header: api-key     E-utilities query: api_key
#   tool + email on every E-utilities call.
#
# JSON is parsed with purrr::pluck() and defaults. NCBI often returns
# length-1 lists or missing fields; never assume a field exists.

.this_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) NULL)
if (is.null(.this_dir) || !file.exists(file.path(.this_dir, "config.R"))) {
  candidates <- c(
    file.path("R", "config.R"),
    "config.R",
    file.path(getwd(), "R", "config.R")
  )
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) stop("Cannot find R/config.R")
  source(hit, local = FALSE)
} else {
  source(file.path(.this_dir, "config.R"), local = FALSE)
}

.ncbi_state <- new.env(parent = emptyenv())
.ncbi_state$api_key <- Sys.getenv("NCBI_API_KEY", unset = "")

#' Set a session-only NCBI API key (sidebar paste). Empty string clears it.
set_session_api_key <- function(key) {
  key <- trimws(as.character(key %||% ""))
  .ncbi_state$api_key <- key
  invisible(key)
}

current_api_key <- function() {
  key <- .ncbi_state$api_key
  if (!nzchar(key)) key <- Sys.getenv("NCBI_API_KEY", unset = "")
  key
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

as_chr1 <- function(x, default = NA_character_) {
  if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) {
    return(default)
  }
  as.character(x[[1]])
}

as_int1 <- function(x, default = NA_integer_) {
  v <- suppressWarnings(as.integer(as_chr1(x, default = NA_character_)))
  if (length(v) == 0 || is.na(v)) default else v
}

ncbi_user_agent <- function() {
  sprintf("%s/1.0 (exploratory teaching app; %s)", NCBI_TOOL, NCBI_EMAIL)
}

throttle_rate <- function() {
  if (nzchar(current_api_key())) 10 else 3
}

is_transient_ncbi <- function(resp) {
  httr2::resp_status(resp) %in% c(429, 500:599)
}

#' Low-level GET; `api_key` is an argument so memoise keys include it.
.ncbi_request <- function(url, query = list(), api_key = current_api_key(),
                          as_text = FALSE) {
  query <- query[!vapply(query, function(z) is.null(z) || (length(z) == 1 && is.na(z)), logical(1))]
  req <- httr2::request(url) |>
    httr2::req_user_agent(ncbi_user_agent()) |>
    httr2::req_url_query(!!!query) |>
    httr2::req_timeout(60) |>
    httr2::req_throttle(rate = throttle_rate()) |>
    httr2::req_retry(
      max_tries = 5,
      retry_on_failure = TRUE,
      is_transient = is_transient_ncbi
    )
  if (nzchar(api_key)) {
    req <- httr2::req_headers(req, `api-key` = api_key)
  }
  resp <- httr2::req_perform(req)
  status <- httr2::resp_status(resp)
  if (status == 404) {
    return(if (as_text) "" else list())
  }
  if (status >= 400) {
    stop(
      "NCBI request failed (", status, ") for ", url,
      ". ", httr2::resp_body_string(resp),
      call. = FALSE
    )
  }
  if (as_text) {
    return(httr2::resp_body_string(resp))
  }
  httr2::resp_body_json(resp, simplifyVector = FALSE)
}

.ncbi_request_memo <- NULL

ncbi_request <- function(url, query = list(), as_text = FALSE) {
  if (is.null(.ncbi_request_memo)) {
    .ncbi_request_memo <<- memoise::memoise(.ncbi_request)
  }
  .ncbi_request_memo(url, query = query, api_key = current_api_key(), as_text = as_text)
}

clear_ncbi_cache <- function() {
  if (!is.null(.ncbi_request_memo)) {
    memoise::forget(.ncbi_request_memo)
  }
  invisible(TRUE)
}

encode_path <- function(...) {
  parts <- unlist(list(...), use.names = FALSE)
  paste(vapply(parts, utils::URLencode, character(1), reserved = TRUE), collapse = "/")
}

datasets_url <- function(...) {
  paste(DATASETS_BASE, encode_path(...), sep = "/")
}

#' Pull every page of a Datasets report that uses page_token / next_page_token.
datasets_paged <- function(url, query = list(), page_size = 1000L) {
  query$page_size <- as.integer(page_size)
  pages <- list()
  token <- NULL
  repeat {
    q <- query
    if (!is.null(token) && nzchar(token)) q$page_token <- token
    json <- ncbi_request(url, query = q)
    pages[[length(pages) + 1L]] <- json
    token <- as_chr1(purrr::pluck(json, "next_page_token"), default = "")
    if (!nzchar(token)) break
    if (length(pages) > 50L) {
      warning("Stopped paging NCBI results after 50 pages.")
      break
    }
  }
  pages
}

# ---- parsers (pure; used by unit tests) ------------------------------------

parse_gene_reports <- function(json) {
  reports <- purrr::pluck(json, "reports", .default = list())
  if (!length(reports)) {
    return(tibble::tibble(
      gene_id = character(), symbol = character(), description = character(),
      tax_id = character(), taxname = character(), common_name = character(),
      type = character(), chromosomes = character(), summary = character(),
      ensembl_gene_ids = character(), swiss_prot = character()
    ))
  }
  purrr::map_dfr(reports, function(item) {
    g <- purrr::pluck(item, "gene", .default = NULL)
    if (is.null(g)) g <- purrr::pluck(item, "product", .default = list())
    chrs <- purrr::pluck(g, "chromosomes", .default = list())
    summary_txt <- purrr::pluck(g, "summary", 1, "description", .default = NA_character_)
    tibble::tibble(
      gene_id = as_chr1(purrr::pluck(g, "gene_id")),
      symbol = as_chr1(purrr::pluck(g, "symbol")),
      description = as_chr1(purrr::pluck(g, "description")),
      tax_id = as_chr1(purrr::pluck(g, "tax_id")),
      taxname = as_chr1(purrr::pluck(g, "taxname")),
      common_name = as_chr1(purrr::pluck(g, "common_name")),
      type = as_chr1(purrr::pluck(g, "type")),
      chromosomes = paste(unlist(chrs), collapse = ","),
      summary = as_chr1(summary_txt),
      ensembl_gene_ids = paste(unlist(purrr::pluck(g, "ensembl_gene_ids", .default = list())), collapse = ","),
      swiss_prot = paste(unlist(purrr::pluck(g, "swiss_prot_accessions", .default = list())), collapse = ",")
    )
  })
}

parse_taxon_suggest <- function(json) {
  hits <- purrr::pluck(json, "sci_name_and_ids", .default = list())
  if (!length(hits)) {
    return(tibble::tibble(
      sci_name = character(), tax_id = character(), common_name = character(),
      matched_term = character(), rank = character()
    ))
  }
  purrr::map_dfr(hits, function(h) {
    tibble::tibble(
      sci_name = as_chr1(purrr::pluck(h, "sci_name")),
      tax_id = as_chr1(purrr::pluck(h, "tax_id")),
      common_name = as_chr1(purrr::pluck(h, "common_name")),
      matched_term = as_chr1(purrr::pluck(h, "matched_term")),
      rank = as_chr1(purrr::pluck(h, "rank"))
    )
  })
}

#' Flatten product_report transcripts into one row per isoform.
parse_product_isoforms <- function(json) {
  reports <- purrr::pluck(json, "reports", .default = list())
  empty <- tibble::tibble(
    gene_id = character(), symbol = character(), tax_id = character(),
    taxname = character(), transcript_acc = character(), transcript_type = character(),
    transcript_length = integer(), cds_length = integer(),
    protein_acc = character(), protein_length = integer(), protein_name = character()
  )
  if (!length(reports)) return(empty)
  rows <- list()
  for (item in reports) {
    prod <- purrr::pluck(item, "product", .default = NULL)
    if (is.null(prod)) prod <- purrr::pluck(item, "gene", .default = list())
    gene_id <- as_chr1(purrr::pluck(prod, "gene_id"))
    symbol <- as_chr1(purrr::pluck(prod, "symbol"))
    tax_id <- as_chr1(purrr::pluck(prod, "tax_id"))
    taxname <- as_chr1(purrr::pluck(prod, "taxname"))
    txs <- purrr::pluck(prod, "transcripts", .default = list())
    if (!length(txs)) {
      rows[[length(rows) + 1L]] <- tibble::tibble(
        gene_id = gene_id, symbol = symbol, tax_id = tax_id, taxname = taxname,
        transcript_acc = NA_character_, transcript_type = NA_character_,
        transcript_length = NA_integer_, cds_length = NA_integer_,
        protein_acc = NA_character_, protein_length = NA_integer_,
        protein_name = NA_character_
      )
      next
    }
    for (tx in txs) {
      cds <- purrr::pluck(tx, "cds", .default = list())
      rng <- purrr::pluck(cds, "range", 1, .default = list())
      cds_len <- NA_integer_
      b <- as_int1(purrr::pluck(rng, "begin"))
      e <- as_int1(purrr::pluck(rng, "end"))
      if (!is.na(b) && !is.na(e)) cds_len <- abs(e - b) + 1L
      prot <- purrr::pluck(tx, "protein", .default = list())
      rows[[length(rows) + 1L]] <- tibble::tibble(
        gene_id = gene_id,
        symbol = symbol,
        tax_id = tax_id,
        taxname = taxname,
        transcript_acc = as_chr1(purrr::pluck(tx, "accession_version")),
        transcript_type = as_chr1(purrr::pluck(tx, "type")),
        transcript_length = as_int1(purrr::pluck(tx, "length")),
        cds_length = cds_len,
        protein_acc = as_chr1(purrr::pluck(prot, "accession_version")),
        protein_length = as_int1(purrr::pluck(prot, "length")),
        protein_name = as_chr1(purrr::pluck(prot, "name"))
      )
    }
  }
  dplyr::bind_rows(rows)
}

#' Prefer NP_/NM_ (curated RefSeq) over XP_/XM_ (models), then longest.
is_curated_refseq <- function(acc, kind = c("protein", "cds")) {
  kind <- match.arg(kind)
  acc <- as.character(acc)
  acc[is.na(acc)] <- ""
  if (kind == "protein") grepl("^(NP_|YP_|AP_|WP_)", acc) else grepl("^(NM_|NR_)", acc)
}

choose_representative_isoform <- function(isoforms, seq_type = c("protein", "cds")) {
  seq_type <- match.arg(seq_type)
  if (is.null(isoforms) || nrow(isoforms) == 0) {
    return(isoforms)
  }
  acc_col <- if (seq_type == "protein") "protein_acc" else "transcript_acc"
  len_col <- if (seq_type == "protein") "protein_length" else "cds_length"
  isoforms |>
    dplyr::filter(!is.na(.data[[acc_col]]), nzchar(.data[[acc_col]])) |>
    dplyr::group_by(gene_id) |>
    dplyr::mutate(
      curated = is_curated_refseq(.data[[acc_col]], kind = seq_type),
      use_len = dplyr::coalesce(.data[[len_col]], transcript_length, 0L)
    ) |>
    dplyr::arrange(dplyr::desc(curated), dplyr::desc(use_len), .data[[acc_col]], .by_group = TRUE) |>
    dplyr::slice(1) |>
    dplyr::ungroup() |>
    dplyr::select(-curated, -use_len)
}

#' Strip .L / .S homeolog suffixes so we can search the same gene in other species.
base_gene_symbol <- function(symbol) {
  s <- as_chr1(symbol, default = "")
  sub("\\.[LS]$", "", s, ignore.case = TRUE)
}

is_xenopus_homeolog <- function(symbol) {
  grepl("\\.[LS]$", as_chr1(symbol, default = ""), ignore.case = TRUE)
}

#' Stable sort then cap: reference genes and X. laevis L/S first, then species name.
cap_orthologs <- function(genes, reference_gene_ids = character(),
                          reference_tax_id = DEFAULT_REFERENCE_TAXID,
                          max_n = MAX_ORTHOLOGS_DEFAULT) {
  max_n <- as.integer(max_n)
  if (max_n > MAX_ORTHOLOGS_HARD) max_n <- MAX_ORTHOLOGS_HARD
  if (is.null(genes) || nrow(genes) == 0) return(genes)
  genes |>
    dplyr::mutate(
      is_ref = gene_id %in% as.character(reference_gene_ids),
      is_xl = tax_id %in% as.character(XENOPUS_LAEVIS_TAXID) |
        grepl("^Xenopus laevis$", taxname, ignore.case = TRUE),
      homeolog_rank = dplyr::case_when(
        grepl("\\.L$", symbol, ignore.case = TRUE) ~ 0L,
        grepl("\\.S$", symbol, ignore.case = TRUE) ~ 1L,
        TRUE ~ 2L
      )
    ) |>
    dplyr::arrange(dplyr::desc(is_ref), dplyr::desc(is_xl), homeolog_rank, taxname, symbol, gene_id) |>
    dplyr::distinct(gene_id, .keep_all = TRUE) |>
    dplyr::slice_head(n = max_n) |>
    dplyr::select(-is_ref, -is_xl, -homeolog_rank)
}

parse_esearch_ids <- function(json) {
  ids <- purrr::pluck(json, "esearchresult", "idlist", .default = list())
  as.character(unlist(ids))
}

parse_fasta_text <- function(text) {
  text <- as_chr1(text, default = "")
  if (!nzchar(trimws(text))) {
    return(tibble::tibble(header = character(), accession = character(), sequence = character()))
  }
  chunks <- strsplit(trimws(text), "\n>")[[1]]
  chunks[1] <- sub("^>", "", chunks[1])
  purrr::map_dfr(chunks, function(chunk) {
    lines <- strsplit(chunk, "\n", fixed = TRUE)[[1]]
    header <- lines[[1]]
    seq <- paste(lines[-1], collapse = "")
    seq <- gsub("\\s+", "", seq)
    acc <- sub("\\s.*$", "", header)
    tibble::tibble(header = header, accession = acc, sequence = seq)
  })
}

# ---- live NCBI calls -------------------------------------------------------

gene_by_symbol <- function(symbol, taxon) {
  symbol <- trimws(as_chr1(symbol, default = ""))
  taxon <- trimws(as_chr1(taxon, default = ""))
  if (!nzchar(symbol) || !nzchar(taxon)) {
    stop("Gene symbol and taxon are required.", call. = FALSE)
  }
  url <- datasets_url("gene", "symbol", symbol, "taxon", taxon)
  pages <- tryCatch(datasets_paged(url), error = function(e) NULL)
  parsed <- if (is.null(pages)) tibble::tibble() else dplyr::bind_rows(lapply(pages, parse_gene_reports))
  if (nrow(parsed) == 0) {
    parsed <- gene_by_symbol_eutils(symbol, taxon)
  }
  parsed
}

gene_by_id <- function(gene_id) {
  gene_id <- trimws(as_chr1(gene_id, default = ""))
  url <- datasets_url("gene", "id", gene_id)
  pages <- tryCatch(datasets_paged(url), error = function(e) NULL)
  parsed <- if (is.null(pages)) tibble::tibble() else dplyr::bind_rows(lapply(pages, parse_gene_reports))
  if (nrow(parsed) == 0) {
    url2 <- datasets_url("gene", "id", gene_id, "dataset_report")
    pages <- tryCatch(datasets_paged(url2), error = function(e) NULL)
    parsed <- if (is.null(pages)) tibble::tibble() else dplyr::bind_rows(lapply(pages, parse_gene_reports))
  }
  parsed
}

fetch_orthologs <- function(gene_id, taxon_filter = NULL) {
  gene_id <- trimws(as_chr1(gene_id, default = ""))
  url <- datasets_url("gene", "id", gene_id, "orthologs")
  query <- list()
  if (!is.null(taxon_filter) && length(taxon_filter)) {
    # Datasets accepts repeated taxon_filter; httr2 expands vectors.
    query$taxon_filter <- as.character(taxon_filter)
  }
  pages <- datasets_paged(url, query = query)
  dplyr::bind_rows(lapply(pages, parse_gene_reports))
}

fetch_product_reports <- function(gene_ids) {
  gene_ids <- unique(as.character(gene_ids))
  gene_ids <- gene_ids[nzchar(gene_ids)]
  if (!length(gene_ids)) {
    return(parse_product_isoforms(list()))
  }
  chunks <- split(gene_ids, ceiling(seq_along(gene_ids) / 10))
  out <- list()
  for (ch in chunks) {
    url <- datasets_url("gene", "id", paste(ch, collapse = ","), "product_report")
    pages <- datasets_paged(url, page_size = 1000L)
    out[[length(out) + 1L]] <- dplyr::bind_rows(lapply(pages, parse_product_isoforms))
  }
  dplyr::bind_rows(out)
}

taxon_suggest <- function(query, restrict_species = TRUE) {
  query <- trimws(as_chr1(query, default = ""))
  if (!nzchar(query) || nchar(query) < 2) {
    return(parse_taxon_suggest(list()))
  }
  url <- datasets_url("taxonomy", "taxon_suggest", query)
  json <- tryCatch(ncbi_request(url, query = list(taxon_resource_filter = "TAXON_RESOURCE_FILTER_GENE")), error = function(e) list())
  tbl <- parse_taxon_suggest(json)
  if (restrict_species && nrow(tbl)) {
    sp <- dplyr::filter(tbl, toupper(rank) == "SPECIES")
    if (nrow(sp)) tbl <- sp
  }
  tbl
}

eutils_query <- function(script, query) {
  api_key <- current_api_key()
  query$tool <- NCBI_TOOL
  query$email <- NCBI_EMAIL
  if (nzchar(api_key)) query$api_key <- api_key
  url <- paste0(EUTILS_BASE, "/", script)
  # E-utilities uses api_key as a query parameter, not the Datasets header.
  req <- httr2::request(url) |>
    httr2::req_user_agent(ncbi_user_agent()) |>
    httr2::req_url_query(!!!query) |>
    httr2::req_timeout(60) |>
    httr2::req_throttle(rate = throttle_rate()) |>
    httr2::req_retry(max_tries = 5, retry_on_failure = TRUE, is_transient = is_transient_ncbi)
  resp <- httr2::req_perform(req)
  if (httr2::resp_status(resp) >= 400) {
    stop("E-utilities failed (", httr2::resp_status(resp), ") for ", script, call. = FALSE)
  }
  if (identical(query$retmode, "json")) {
    httr2::resp_body_json(resp, simplifyVector = FALSE)
  } else {
    httr2::resp_body_string(resp)
  }
}

gene_by_symbol_eutils <- function(symbol, taxon) {
  term <- sprintf("%s[sym] AND %s[taxid]", symbol, taxon)
  json <- tryCatch(
    eutils_query("esearch.fcgi", list(db = "gene", term = term, retmode = "json", retmax = 20)),
    error = function(e) NULL
  )
  ids <- if (is.null(json)) character() else parse_esearch_ids(json)
  if (!length(ids)) {
    return(parse_gene_reports(list()))
  }
  gene_by_id(paste(ids, collapse = ","))
}

chunk_vec <- function(x, n = 50L) {
  if (!length(x)) return(list())
  split(x, ceiling(seq_along(x) / n))
}

#' Download FASTA. Never mix protein and nuccore accessions in one call.
efetch_fasta <- function(accessions, db = c("protein", "nuccore")) {
  db <- match.arg(db)
  accessions <- unique(as.character(accessions))
  accessions <- accessions[nzchar(accessions) & !is.na(accessions)]
  if (!length(accessions)) {
    return(parse_fasta_text(""))
  }
  parts <- list()
  for (ch in chunk_vec(accessions, 50L)) {
    txt <- eutils_query(
      "efetch.fcgi",
      list(db = db, id = paste(ch, collapse = ","), rettype = "fasta", retmode = "text")
    )
    parts[[length(parts) + 1L]] <- parse_fasta_text(txt)
  }
  dplyr::bind_rows(parts)
}

#' Resolve a typed gene (symbol or GeneID) in the reference taxon.
resolve_reference_genes <- function(query, taxon) {
  query <- trimws(as_chr1(query, default = ""))
  if (!nzchar(query)) stop("Enter a gene symbol or NCBI Gene ID.", call. = FALSE)
  if (grepl("^[0-9]+$", query)) {
    genes <- gene_by_id(query)
  } else {
    genes <- gene_by_symbol(query, taxon)
  }
  genes <- genes[!is.na(genes$gene_id) & nzchar(genes$gene_id), , drop = FALSE]
  if (!nrow(genes)) {
    stop(
      "Gene not found for '", query, "' in taxon ", taxon,
      ". Check the spelling or try a numeric NCBI Gene ID.",
      call. = FALSE
    )
  }
  genes
}

#' Build the working ortholog table: NCBI orthologs, plus homeologs, plus
#' per-species symbol search when the ortholog endpoint is empty (common for
#' some Xenopus genes).
collect_ortholog_table <- function(reference_genes, comparison_taxids = NULL,
                                   use_all_orthologs = FALSE,
                                   max_n = MAX_ORTHOLOGS_DEFAULT) {
  ref_ids <- unique(as.character(reference_genes$gene_id))
  ref_tax <- as_chr1(reference_genes$tax_id[1], default = DEFAULT_REFERENCE_TAXID)
  symbol0 <- base_gene_symbol(reference_genes$symbol[1])

  ortho <- tibble::tibble()
  for (gid in ref_ids) {
    got <- tryCatch(fetch_orthologs(gid), error = function(e) tibble::tibble())
    ortho <- dplyr::bind_rows(ortho, got)
  }

  # Some Xenopus genes are not in a large NCBI Ortholog set. Fall back to the
  # human (or first comparison) gene of the same symbol, then merge.
  if (nrow(ortho) < 3) {
    fallback_tax <- c("9606", as.character(comparison_taxids))
    fallback_tax <- unique(fallback_tax[nzchar(fallback_tax)])
    for (tax in fallback_tax) {
      extra <- tryCatch(gene_by_symbol(symbol0, tax), error = function(e) tibble::tibble())
      if (nrow(extra)) {
        hid <- extra$gene_id[[1]]
        more <- tryCatch(fetch_orthologs(hid), error = function(e) extra)
        ortho <- dplyr::bind_rows(ortho, extra, more)
        if (nrow(dplyr::distinct(ortho, gene_id)) >= 3) break
      }
    }
  }

  # Always keep the reference genes (both X. laevis homeologs when present).
  ortho <- dplyr::bind_rows(reference_genes, ortho) |>
    dplyr::filter(!is.na(gene_id), nzchar(gene_id)) |>
    dplyr::distinct(gene_id, .keep_all = TRUE)

  if (!use_all_orthologs && length(comparison_taxids)) {
    keep_tax <- unique(c(ref_tax, as.character(comparison_taxids)))
    # Symbol search for any selected species still missing.
    present <- unique(ortho$tax_id)
    missing_tax <- setdiff(keep_tax, present)
    for (tax in missing_tax) {
      extra <- tryCatch(gene_by_symbol(symbol0, tax), error = function(e) tibble::tibble())
      ortho <- dplyr::bind_rows(ortho, extra)
    }
    ortho <- dplyr::bind_rows(reference_genes, ortho) |>
      dplyr::filter(tax_id %in% keep_tax | gene_id %in% ref_ids) |>
      dplyr::distinct(gene_id, .keep_all = TRUE)
  }

  if (!nrow(ortho)) {
    stop("No orthologs were returned by NCBI for this gene.", call. = FALSE)
  }

  cap_orthologs(
    ortho,
    reference_gene_ids = ref_ids,
    reference_tax_id = ref_tax,
    max_n = max_n
  )
}
