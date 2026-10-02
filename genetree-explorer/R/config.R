# Lab-editable NCBI identity for E-utilities.
# NCBI asks callers to send a tool name and a contact email. Change these
# two values to match your group. They are sent only to NCBI, never shown
# in the Shiny user interface.

NCBI_TOOL <- "GeneTreeExplorer"
NCBI_EMAIL <- "mcfarlane-lab@ucalgary.ca"

DATASETS_BASE <- "https://api.ncbi.nlm.nih.gov/datasets/v2"
EUTILS_BASE <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Hard cap on how many ortholog tips we will align. The UI slider cannot
# go above this value (spec: default 12, hard max 40).
MAX_ORTHOLOGS_HARD <- 40L
MAX_ORTHOLOGS_DEFAULT <- 12L
ML_MAX_SEQUENCES <- 12L
BOOTSTRAP_CONFIRM_N <- 20L

# Default comparison species (NCBI Taxonomy IDs).
DEFAULT_SPECIES <- c(
  "Homo sapiens" = "9606",
  "Mus musculus" = "10090",
  "Gallus gallus" = "9031",
  "Danio rerio" = "7955",
  "Xenopus tropicalis" = "8364",
  "Xenopus laevis" = "8355"
)

DEFAULT_REFERENCE_TAXID <- "8355"
DEFAULT_REFERENCE_NAME <- "Xenopus laevis"
XENOPUS_LAEVIS_TAXID <- "8355"
