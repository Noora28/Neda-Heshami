# GeneTreeExplorer

Exploratory Shiny app: type a gene symbol, pick species, download NCBI
orthologs, align **one representative isoform per gene**, and draw a tree.

This is a **teaching / exploration** tool. It is **not** Ensembl Compara,
OrthoFinder, or a publication phylogeny. Trees are a species-sampled gene
tree from NCBI Gene orthologs. Do not publish them without further analysis.

The app lives in `genetree-explorer/`. The repository root still holds the
lab profile `index.html` / `style.css`.

## NCBI endpoints used (checked 2 Oct 2026)

**NCBI Datasets REST API v2** base: `https://api.ncbi.nlm.nih.gov/datasets/v2`

Documented in the NCBI Datasets OpenAPI spec
(`https://github.com/ncbi/datasets/blob/master/datasets.openapi.yaml`) and
`https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/`.

| Purpose | Method | Path | JSON notes |
| --- | --- | --- | --- |
| Gene by symbol + taxon | GET | `/gene/symbol/{symbols}/taxon/{taxon}` | `reports[].gene` |
| Gene by ID | GET | `/gene/id/{gene_ids}` | Works; OpenAPI marks it deprecated. Current alias: `/gene/id/{gene_ids}/dataset_report` |
| Orthologs | GET | `/gene/id/{gene_id}/orthologs` | `reports[].gene`; `page_token` / `next_page_token` |
| Products / isoforms | GET | `/gene/id/{gene_ids}/product_report` | **`reports[].product`** (not `.gene`); isoforms in `product.transcripts[]` with `protein.accession_version` |
| Taxonomy autocomplete | GET | `/taxonomy/taxon_suggest/{taxon_query}` | `sci_name_and_ids[]` (`sci_name`, `tax_id`, `common_name`, `rank`) |

**E-utilities** base: `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/`
(Book: https://www.ncbi.nlm.nih.gov/books/NBK25499/ , last update 4 Mar 2026)

- `esearch.fcgi` / `esummary.fcgi` with `retmode=json` (gene fallback)
- `efetch.fcgi` `db=protein` or `db=nuccore`, `rettype=fasta`, `retmode=text`
- Datasets header: `api-key`. E-utilities query: `api_key`. Always send `tool` and `email`.

Edit the placeholder lab email in `R/config.R` (`NCBI_EMAIL`). It is **not**
shown in the UI.

Some Xenopus genes (including `pax6.L`) return only themselves from the
orthologs endpoint. The app then searches the same base symbol in the
comparison species and, if needed, uses the human gene’s ortholog set.

## Package pin

- CRAN mirror: `https://cloud.r-project.org`
- Bioconductor: **3.18** on R 4.3, **3.19** on R 4.4, **3.21** on R 4.5+
  (`install_packages.R` chooses this automatically).
- No IQ-TREE or other external binaries.

A generated `renv.lock` is included when `renv::snapshot()` succeeds after
install. If it is missing, `install_packages.R` plus this pin is enough.

## Windows: run from Neda’s machine

Target folder:

`C:\Neda's work\Webapp\phylotree\genetree-explorer`

1. Install [R](https://cran.r-project.org/) 4.3 or newer (64-bit).
2. Optional but useful on Windows: [RTools](https://cran.r-project.org/bin/windows/Rtools/)
   only if Bioconductor has to compile packages from source. Binary packages
   are used when available.
3. Open **R** or **RStudio**.
4. Install packages once (takes several minutes):

```r
setwd("C:/Neda's work/Webapp/phylotree/genetree-explorer")
source("install_packages.R")
```

5. Optional NCBI API key (never commit it). In that same folder create
   `.Renviron` (Notepad) with one line:

```
NCBI_API_KEY=paste_your_key_here
```

   You can also paste a key in the sidebar for the current session only.

6. Start the app:

```r
setwd("C:/Neda's work/Webapp/phylotree/genetree-explorer")
shiny::runApp()
```

   In RStudio: File → Open Folder → that path → open `app.R` → Run App.

7. In the browser: gene `pax6`, reference *Xenopus laevis*, default comparison
   species, protein, NJ, bootstrap 0 → **Run analysis**.

PowerShell equivalent:

```powershell
Set-Location "C:\Neda's work\Webapp\phylotree\genetree-explorer"
Rscript -e "shiny::runApp(getwd())"
```

## Tests

```r
setwd("C:/Neda's work/Webapp/phylotree/genetree-explorer")
source("install_packages.R")   # once
Rscript tests/run_tests.R
```

Optional live NCBI smoke (skips / records failure; does not block):

```r
Rscript tests/live_smoke.R
```

### Live NCBI smoke result

Recorded during development (2 Oct 2026): see the “Live NCBI smoke” section
at the bottom of this file. Re-run `tests/live_smoke.R` on the lab machine
to refresh.

## v1 analysis path

Protein sequences, neighbour-joining, midpoint root (unless you name a unique
outgroup tip). That path is the default and should work whenever NCBI returns
at least three sequences. UPGMA and ML are optional; ML is refused for n > 12
and never bootstraps.
