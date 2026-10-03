# GeneTreeExplorer Shiny app
# Run from this folder:  shiny::runApp()
# Fetch from NCBI happens only when the user clicks Run.

app_dir <- getwd()
source(file.path(app_dir, "R", "ncbi_api.R"), local = FALSE)
source(file.path(app_dir, "R", "phylo.R"), local = FALSE)

disclaimer_text <- paste(
  "Teaching / exploration tool: this is a species-sampled gene tree built from",
  "NCBI Gene orthologs and one representative isoform per gene. It is not",
  "Ensembl Compara, OrthoFinder, or a publication phylogeny. Do not publish",
  "these trees without further analysis."
)

default_species_choices <- DEFAULT_SPECIES

ui <- bslib::page_sidebar(
  title = "GeneTreeExplorer",
  theme = bslib::bs_theme(version = 5, bootswatch = "flatly"),
  sidebar = bslib::sidebar(
    width = 360,
    shiny::p(disclaimer_text, class = "small text-muted"),
    shiny::textInput(
      "gene_query",
      "Gene symbol or NCBI Gene ID",
      value = "pax6",
      placeholder = "e.g. pax6 or 379100"
    ),
    shiny::selectizeInput(
      "ref_species",
      "Reference species",
      choices = default_species_choices,
      selected = DEFAULT_REFERENCE_TAXID,
      options = list(create = FALSE, placeholder = "Type a species name")
    ),
    shiny::selectizeInput(
      "comp_species",
      "Comparison species",
      choices = default_species_choices,
      selected = unname(default_species_choices),
      multiple = TRUE,
      options = list(placeholder = "Select one or more species")
    ),
    shiny::textInput(
      "tax_search",
      "Search NCBI taxonomy to add a species",
      placeholder = "e.g. zebrafish"
    ),
    shiny::helpText("Autocomplete uses Datasets GET /taxonomy/taxon_suggest/{query}."),
    shiny::checkboxInput(
      "use_all_orthologs",
      "Use all orthologs returned by NCBI",
      value = FALSE
    ),
    shiny::sliderInput(
      "max_orthologs",
      "Maximum orthologs (hard max 40)",
      min = 3, max = MAX_ORTHOLOGS_HARD, value = MAX_ORTHOLOGS_DEFAULT, step = 1
    ),
    shiny::radioButtons(
      "seq_type",
      "Sequence type",
      choices = c("Protein (recommended)" = "protein", "CDS (coding DNA)" = "cds"),
      selected = "protein"
    ),
    shiny::conditionalPanel(
      "input.seq_type == 'cds'",
      shiny::div(
        class = "alert alert-warning small",
        "CDS is aligned by translating to protein, aligning amino acids, then back-translating. Sequences without a complete coding frame (length not divisible by 3) are refused."
      )
    ),
    shiny::radioButtons(
      "tree_method",
      "Tree method",
      choices = c(
        "Neighbour-joining (default)" = "NJ",
        "UPGMA" = "UPGMA",
        "Maximum likelihood (no bootstrap)" = "ML"
      ),
      selected = "NJ"
    ),
    shiny::conditionalPanel(
      "input.tree_method == 'UPGMA'",
      shiny::div(
        class = "alert alert-warning small",
        "UPGMA assumes a molecular clock (the same substitution rate on every branch). That is often unrealistic for these genes."
      )
    ),
    shiny::selectInput(
      "bootstrap",
      "Bootstrap replicates",
      choices = c("0 (none)" = 0, "100" = 100),
      selected = 0
    ),
    shiny::selectizeInput(
      "outgroup",
      "Outgroup tip (optional)",
      choices = c("None (midpoint root)" = ""),
      selected = ""
    ),
    shiny::helpText(
      "Default is midpoint root (no outgroup). Do not root on human when the question is about vertebrate ancestry. A distant species such as zebrafish is a better optional outgroup, but it is not required."
    ),
    shiny::hr(),
    shiny::passwordInput(
      "api_key",
      "NCBI API key (session only, optional)",
      placeholder = "Not stored; also reads NCBI_API_KEY from .Renviron"
    ),
    shiny::helpText("3 requests/second without a key, 10/second with one. Never commit a key to git."),
    shiny::actionButton("run", "Run analysis", class = "btn-primary w-100"),
    shiny::actionButton("rebuild", "Rebuild tree from selected orthologs", class = "btn-outline-secondary w-100")
  ),
  bslib::navset_card_tab(
    id = "main_tabs",
    bslib::nav_panel("Gene information", DT::DTOutput("gene_table") |> shinycssloaders::withSpinner()),
    bslib::nav_panel("Orthologs", shiny::div(
      shiny::p("Select rows to include in the alignment and tree. Both Xenopus laevis L and S homeologs are kept as separate tips when NCBI returns them."),
      DT::DTOutput("ortholog_table") |> shinycssloaders::withSpinner()
    )),
    bslib::nav_panel("Sequences", shiny::verbatimTextOutput("seq_fasta") |> shinycssloaders::withSpinner()),
    bslib::nav_panel("Alignment", shiny::div(
      shiny::textOutput("aln_identity"),
      shiny::verbatimTextOutput("aln_fasta")
    )),
    bslib::nav_panel("Tree", shiny::div(
      shiny::radioButtons(
        "tree_layout",
        "Layout",
        choices = c("Rectangular" = "rectangular", "Circular" = "circular"),
        selected = "rectangular",
        inline = TRUE
      ),
      shiny::checkboxInput("show_bootstrap", "Show bootstrap values", value = TRUE),
      shiny::helpText("Bootstrap values are hidden on circular layouts (they overlap tips)."),
      shiny::plotOutput("tree_plot", height = "640px") |> shinycssloaders::withSpinner(),
      shiny::verbatimTextOutput("tree_status")
    )),
    bslib::nav_panel("How to read this tree", shiny::div(
      shiny::h4("What this tree is"),
      shiny::p("Each tip is one gene from one species (two tips for Xenopus laevis L and S homeologs). Branches show similarity of the aligned representative isoforms, not speciation times."),
      shiny::h4("Rooting"),
      shiny::p("With no outgroup, the tree is midpoint-rooted (the root is placed halfway along the longest tip-to-tip path). That is a convenience, not biological evidence of the ancestor."),
      shiny::h4("What to distrust"),
      shiny::tags$ul(
        shiny::tags$li("One isoform per gene can miss biologically important splice forms."),
        shiny::tags$li("NCBI orthologs are not always complete, especially in amphibians."),
        shiny::tags$li("NJ is fast and exploratory. UPGMA assumes a clock. ML without bootstrap has no support values."),
        shiny::tags$li("Do not treat this figure as a publication phylogeny.")
      )
    )),
    bslib::nav_panel("Raw JSON", shiny::verbatimTextOutput("raw_json")),
    bslib::nav_panel("Downloads", shiny::div(
      shiny::downloadButton("dl_gene", "Gene CSV"),
      shiny::downloadButton("dl_ortho", "Ortholog CSV"),
      shiny::downloadButton("dl_seq", "Sequences FASTA"),
      shiny::downloadButton("dl_aln", "Alignment FASTA"),
      shiny::downloadButton("dl_nwk", "Tree Newick"),
      shiny::downloadButton("dl_png", "Tree PNG"),
      shiny::downloadButton("dl_pdf", "Tree PDF")
    ))
  )
)

server <- function(input, output, session) {
  rv <- shiny::reactiveValues(
    gene = NULL,
    orthologs = NULL,
    isoforms = NULL,
    seqs = NULL,
    aln = NULL,
    tree = NULL,
    tree_msg = NULL,
    identity = NA_real_,
    raw_json = "",
    status = "Idle. Click Run analysis to fetch from NCBI."
  )

  # Taxonomy autocomplete for reference + comparison species.
  rv$species_choices <- default_species_choices

  shiny::observeEvent(input$tax_search, {
    q <- trimws(input$tax_search)
    if (nchar(q) < 3) return()
    hits <- tryCatch(taxon_suggest(q), error = function(e) NULL)
    if (is.null(hits) || !nrow(hits)) return()
    extra <- stats::setNames(
      hits$tax_id,
      paste0(hits$sci_name, " (", hits$common_name, ")")
    )
    rv$species_choices <- c(rv$species_choices, extra)
    rv$species_choices <- rv$species_choices[!duplicated(unname(rv$species_choices))]
    shiny::updateSelectizeInput(
      session, "ref_species",
      choices = rv$species_choices,
      selected = input$ref_species,
      server = TRUE
    )
    shiny::updateSelectizeInput(
      session, "comp_species",
      choices = rv$species_choices,
      selected = input$comp_species,
      server = TRUE
    )
  }, ignoreInit = TRUE)

  shiny::observeEvent(input$tree_method, {
    method <- input$tree_method
    if (identical(method, "ML")) {
      shiny::updateSelectInput(
        session, "bootstrap",
        choices = c("0 (ML never bootstraps)" = 0),
        selected = 0
      )
    } else {
      shiny::updateSelectInput(
        session, "bootstrap",
        choices = c("0 (none)" = 0, "100" = 100, "500" = 500, "1000" = 1000),
        selected = 0
      )
    }
  })

  shiny::observeEvent(input$api_key, {
    set_session_api_key(input$api_key)
  }, ignoreInit = TRUE)

  selected_orthologs <- shiny::reactive({
    shiny::req(rv$orthologs)
    sel <- input$ortholog_table_rows_selected
    if (is.null(sel) || !length(sel)) rv$orthologs else rv$orthologs[sel, , drop = FALSE]
  })

  run_pipeline <- function(orthologs_already = NULL) {
    set_session_api_key(input$api_key)
    gene_q <- trimws(input$gene_query)
    if (!nzchar(gene_q)) stop("Enter a gene symbol or NCBI Gene ID.", call. = FALSE)
    ref_tax <- as.character(input$ref_species)
    if (!nzchar(ref_tax)) stop("Choose a reference species.", call. = FALSE)

    shiny::withProgress(message = "GeneTreeExplorer", value = 0, {
      incProgress <- function(amount, detail) shiny::incProgress(amount = amount, detail = detail)

      incProgress(0.05, "Resolving gene at NCBI")
      gene_tbl <- resolve_reference_genes(gene_q, ref_tax)
      rv$gene <- gene_tbl
      rv$raw_json <- jsonlite::toJSON(as.list(gene_tbl), pretty = TRUE, auto_unbox = TRUE)

      if (is.null(orthologs_already)) {
        incProgress(0.15, "Fetching orthologs")
        use_all <- isTRUE(input$use_all_orthologs)
        max_n <- as.integer(input$max_orthologs)
        comp <- if (use_all) NULL else as.character(input$comp_species)
        ortho <- collect_ortholog_table(
          gene_tbl,
          comparison_taxids = comp,
          use_all_orthologs = use_all,
          max_n = max_n
        )
        rv$orthologs <- ortho
      } else {
        ortho <- orthologs_already
        incProgress(0.15, "Using selected orthologs")
      }

      n_ortho <- nrow(ortho)
      if (!n_ortho) stop("No orthologs were returned by NCBI.", call. = FALSE)

      incProgress(0.15, "Choosing one isoform per gene")
      isoforms <- fetch_product_reports(ortho$gene_id)
      chosen <- choose_representative_isoform(isoforms, seq_type = input$seq_type)
      chosen <- dplyr::left_join(
        ortho,
        chosen[, intersect(names(chosen), c(
          "gene_id", "transcript_acc", "transcript_type", "transcript_length",
          "cds_length", "protein_acc", "protein_length", "protein_name"
        )), drop = FALSE],
        by = "gene_id"
      )
      rv$isoforms <- chosen

      incProgress(0.2, "Downloading sequences (E-utilities efetch)")
      seq_type <- input$seq_type
      if (identical(seq_type, "protein")) {
        accs <- chosen$protein_acc
        fasta <- efetch_fasta(accs, db = "protein")
      } else {
        accs <- chosen$transcript_acc
        fasta <- efetch_fasta(accs, db = "nuccore")
      }
      seqs <- sequences_from_fasta_tbl(fasta, chosen, seq_type)
      seqs <- seqs[!is.na(seqs$sequence) & nzchar(seqs$sequence), , drop = FALSE]
      if (!nrow(seqs)) {
        stop("Invalid sequences: NCBI FASTA download returned nothing we could parse.", call. = FALSE)
      }
      rv$seqs <- seqs
      shiny::updateSelectizeInput(
        session, "outgroup",
        choices = c("None (midpoint root)" = "", setNames(seqs$tip_label, seqs$tip_label)),
        selected = ""
      )

      n <- nrow(seqs)
      few <- skip_tree_reason(n)
      if (!is.null(few)) {
        rv$aln <- NULL
        rv$tree <- NULL
        rv$tree_msg <- few
        rv$identity <- NA_real_
        rv$status <- few
        return(invisible(NULL))
      }

      method <- input$tree_method
      if (identical(method, "ML")) {
        mlr <- ml_refused_reason(n)
        if (!is.null(mlr)) {
          rv$aln <- NULL
          rv$tree <- NULL
          rv$tree_msg <- mlr
          rv$status <- mlr
          return(invisible(NULL))
        }
      }

      bs <- as.integer(input$bootstrap)
      if (identical(method, "ML")) bs <- 0L
      if (bs > 0 && n > BOOTSTRAP_CONFIRM_N) {
        rv$status <- paste0(
          "Bootstrap with n > ", BOOTSTRAP_CONFIRM_N,
          " can take several minutes. Still working…"
        )
      }

      incProgress(0.2, "Aligning sequences")
      aln <- align_sequences(seqs, seq_type = seq_type)
      rv$aln <- aln
      rv$identity <- mean_pairwise_identity(aln)

      incProgress(0.2, if (bs > 0) "Building tree (bootstrap still working…)" else "Building tree")
      built <- build_tree(
        aln,
        seq_type = seq_type,
        method = method,
        bootstrap = bs,
        outgroup_tip = input$outgroup
      )
      rv$tree <- built$tree
      rv$tree_msg <- built$message
      rv$status <- if (is.null(built$message)) {
        paste0("Done. ", n, " sequences, method ", method, ", bootstrap ", bs, ".")
      } else built$message
    })
  }

  shiny::observeEvent(input$run, {
    tryCatch(
      run_pipeline(NULL),
      error = function(e) {
        msg <- conditionMessage(e)
        if (grepl("429|rate", msg, ignore.case = TRUE)) {
          msg <- paste("NCBI rate limit:", msg, "Wait a few seconds or add an API key.")
        } else if (grepl("Failed to fetch|Could not resolve|timed out|timeout|502|503", msg, ignore.case = TRUE)) {
          msg <- paste("NCBI looks down or unreachable:", msg)
        }
        rv$status <- msg
        rv$tree_msg <- msg
        shiny::showNotification(msg, type = "error", duration = 12)
      }
    )
  })

  shiny::observeEvent(input$rebuild, {
    shiny::req(rv$orthologs)
    tryCatch(
      run_pipeline(selected_orthologs()),
      error = function(e) {
        rv$status <- conditionMessage(e)
        shiny::showNotification(conditionMessage(e), type = "error", duration = 12)
      }
    )
  })

  output$gene_table <- DT::renderDT({
    shiny::validate(shiny::need(!is.null(rv$gene), rv$status %||% "Click Run analysis."))
    DT::datatable(rv$gene, rownames = FALSE, options = list(scrollX = TRUE))
  })

  output$ortholog_table <- DT::renderDT({
    shiny::validate(shiny::need(!is.null(rv$orthologs), "Click Run analysis."))
    DT::datatable(
      rv$orthologs,
      rownames = FALSE,
      selection = list(mode = "multiple", selected = seq_len(nrow(rv$orthologs))),
      options = list(scrollX = TRUE, pageLength = 25)
    )
  })

  output$seq_fasta <- shiny::renderText({
    shiny::validate(shiny::need(!is.null(rv$seqs), "Click Run analysis."))
    ss <- rv$seqs$sequence
    names(ss) <- rv$seqs$tip_label
    unaligned_to_fasta_text(ss)
  })

  output$aln_identity <- shiny::renderText({
    shiny::validate(shiny::need(!is.null(rv$aln), "Click Run analysis."))
    sprintf(
      "Mean pairwise identity over ungapped positions: %.1f%%",
      100 * rv$identity
    )
  })

  output$aln_fasta <- shiny::renderText({
    shiny::validate(shiny::need(!is.null(rv$aln), "Click Run analysis."))
    alignment_to_fasta_text(rv$aln)
  })

  output$tree_status <- shiny::renderText({
    paste(rv$status, rv$tree_msg, sep = "\n")
  })

  output$tree_plot <- shiny::renderPlot({
    shiny::validate(shiny::need(!is.null(rv$tree) || !is.null(rv$tree_msg), "Click Run analysis."))
    if (is.null(rv$tree)) {
      plot.new()
      title(main = rv$tree_msg %||% "No tree")
      return(invisible(NULL))
    }
    show_bs <- isTRUE(input$show_bootstrap) && !identical(input$tree_layout, "circular")
    print(plot_ggtree(rv$tree, layout = input$tree_layout, show_bootstrap = show_bs))
  })

  output$raw_json <- shiny::renderText({
    rv$raw_json %||% ""
  })

  output$dl_gene <- shiny::downloadHandler(
    filename = function() "genetree_gene.csv",
    content = function(file) utils::write.csv(rv$gene, file, row.names = FALSE)
  )
  output$dl_ortho <- shiny::downloadHandler(
    filename = function() "genetree_orthologs.csv",
    content = function(file) utils::write.csv(rv$orthologs, file, row.names = FALSE)
  )
  output$dl_seq <- shiny::downloadHandler(
    filename = function() "genetree_sequences.fasta",
    content = function(file) writeLines(unaligned_to_fasta_text(setNames(rv$seqs$sequence, rv$seqs$tip_label)), file)
  )
  output$dl_aln <- shiny::downloadHandler(
    filename = function() "genetree_alignment.fasta",
    content = function(file) writeLines(alignment_to_fasta_text(rv$aln), file)
  )
  output$dl_nwk <- shiny::downloadHandler(
    filename = function() "genetree_tree.nwk",
    content = function(file) ape::write.tree(rv$tree, file = file)
  )
  output$dl_png <- shiny::downloadHandler(
    filename = function() "genetree_tree.png",
    content = function(file) {
      grDevices::png(file, width = 1400, height = 1000, res = 120)
      print(plot_ggtree(rv$tree, layout = input$tree_layout, show_bootstrap = isTRUE(input$show_bootstrap)))
      grDevices::dev.off()
    }
  )
  output$dl_pdf <- shiny::downloadHandler(
    filename = function() "genetree_tree.pdf",
    content = function(file) {
      grDevices::pdf(file, width = 11, height = 8.5)
      print(plot_ggtree(rv$tree, layout = input$tree_layout, show_bootstrap = isTRUE(input$show_bootstrap)))
      grDevices::dev.off()
    }
  )
}

shiny::shinyApp(ui, server)
