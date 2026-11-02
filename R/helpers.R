# helpers.R -----------------------------------------------------------------
# Input parsing + small UI builders used by app.R.
# ---------------------------------------------------------------------------

# Parse free-text input into a data.frame(SMILES, ID).
# Each non-empty line: "<SMILES>" or "<SMILES><sep><ID>" where sep is a comma,
# tab, or run of spaces. Lines starting with '#' are ignored.
parse_smiles_text <- function(txt) {
  if (is.null(txt) || !nzchar(trimws(txt))) return(empty_compounds())
  lines <- strsplit(txt, "\r?\n")[[1]]
  lines <- trimws(lines)
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  if (!length(lines)) return(empty_compounds())

  smi <- character(0); id <- character(0)
  for (ln in lines) {
    parts <- strsplit(ln, "[,\t]|\\s{1,}", perl = TRUE)[[1]]
    parts <- parts[nzchar(parts)]
    smi <- c(smi, parts[1])
    id  <- c(id, if (length(parts) >= 2) paste(parts[-1], collapse = "_") else NA_character_)
  }
  finalize_compounds(smi, id)
}

# Parse an uploaded file (.smi/.txt/.csv) into compounds.
parse_smiles_file <- function(path, name = "") {
  ext <- tolower(tools::file_ext(name))
  if (ext == "csv") {
    df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    cn <- tolower(names(df))
    scol <- which(grepl("smiles|structure", cn))[1]
    if (is.na(scol)) scol <- 1
    icol <- which(grepl("^id$|name|molecule", cn))[1]
    smi <- as.character(df[[scol]])
    id  <- if (!is.na(icol)) as.character(df[[icol]]) else NA_character_
    return(finalize_compounds(smi, id))
  }
  # .smi / .txt : one record per line, SMILES then optional id
  parse_smiles_text(paste(readLines(path, warn = FALSE), collapse = "\n"))
}

empty_compounds <- function() {
  data.frame(SMILES = character(0), ID = character(0), stringsAsFactors = FALSE)
}

finalize_compounds <- function(smi, id) {
  smi <- trimws(as.character(smi))
  keep <- nzchar(smi)
  smi <- smi[keep]; id <- id[keep]
  if (!length(smi)) return(empty_compounds())
  id <- as.character(id)
  blank <- is.na(id) | !nzchar(trimws(id))
  id[blank] <- paste0("Mol_", which(blank))
  id <- make.unique(id, sep = "_")
  data.frame(SMILES = smi, ID = id, stringsAsFactors = FALSE)
}

# --- small presentation helpers -------------------------------------------

# Applicability-domain badge (in / out / unknown).
ad_badge <- function(ad) {
  if (is.na(ad)) return(htmltools::span(class = "badge bg-secondary", "AD: n/a"))
  if (ad >= 0.5) htmltools::span(class = "badge bg-success", "In domain")
  else           htmltools::span(class = "badge bg-warning text-dark", "Out of domain")
}

# Confidence index -> bootstrap colour.
conf_class <- function(conf) {
  if (is.na(conf)) return("bg-secondary")
  if (conf >= 0.7) "bg-success" else if (conf >= 0.5) "bg-info" else "bg-warning text-dark"
}

fmt_num <- function(x, digits = 3) {
  if (is.na(x)) return("—")
  formatC(x, format = "fg", digits = digits, big.mark = ",")
}
