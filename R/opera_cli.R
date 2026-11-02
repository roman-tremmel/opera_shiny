# opera_cli.R ---------------------------------------------------------------
# Thin wrapper that shells out to a locally installed OPERA command-line
# executable, runs the three ADME models, and returns the parsed predictions.
#
# OPERA CLI reference (v2.6-2.9). Core flags used here:
#   -i  <file>        input structure file (.smi / .sdf / .mol / .csv / .txt)
#   -o  <file>        output file (.csv)
#   -e  <list>        endpoints to predict (e.g. CACO2 FUB Clint)
#   -v  <0|1|2>       verbosity
#   -n                append nearest-neighbour info      (optional)
#   -x                append experimental values if known (optional)
#   -c                clean temporary files               (optional)
# The exact executable and flags are configurable so the app matches the
# user's installed build.
# ---------------------------------------------------------------------------

# Resolve the OPERA command: explicit arg > env var OPERA_CMD > "OPERA".
opera_command <- function(cmd = NULL) {
  if (!is.null(cmd) && nzchar(cmd)) return(cmd)
  env <- Sys.getenv("OPERA_CMD", unset = "")
  if (nzchar(env)) return(env)
  "OPERA"
}

# Is an OPERA executable actually reachable?
opera_available <- function(cmd = NULL) {
  cmd <- opera_command(cmd)
  # Absolute/relative path that exists, or something on PATH.
  if (file.exists(cmd)) return(TRUE)
  nzchar(Sys.which(cmd))
}

# Write a data.frame of compounds (SMILES + ID) to a .smi file OPERA can read.
# OPERA .smi format: "<SMILES>\t<ID>" per line, no header.
write_smi <- function(compounds, path) {
  df <- compounds
  ids <- if (!is.null(df$ID)) df$ID else paste0("Mol_", seq_len(nrow(df)))
  lines <- paste(df$SMILES, ids, sep = "\t")
  writeLines(lines, path)
  invisible(path)
}

#' Run OPERA's ADME models on a set of compounds.
#'
#' @param compounds data.frame with columns SMILES (chr) and optional ID (chr)
#' @param cmd       path to OPERA executable (NULL -> resolve via env/PATH)
#' @param endpoints OPERA endpoint tokens to request
#' @param neighbors append nearest neighbours (-n)
#' @param exp       append experimental values (-x)
#' @param verbose   verbosity level for -v
#' @param timeout   seconds before the call is aborted
#' @return list(ok, tidy, raw, command, log, error)
run_opera <- function(compounds,
                      cmd       = NULL,
                      endpoints = adme_tokens(),
                      neighbors = FALSE,
                      exp       = FALSE,
                      verbose   = 1L,
                      timeout   = 1800) {

  cmd <- opera_command(cmd)
  res <- list(ok = FALSE, tidy = NULL, raw = NULL,
              command = NA_character_, log = character(0), error = NULL)

  if (!nrow(compounds)) {
    res$error <- "No compounds supplied."
    return(res)
  }
  if (!opera_available(cmd)) {
    res$error <- paste0(
      "OPERA executable not found (tried '", cmd, "'). Install OPERA and set ",
      "its path in the sidebar, or the OPERA_CMD environment variable."
    )
    return(res)
  }

  workdir <- tempfile("opera_run_")
  dir.create(workdir)
  on.exit(unlink(workdir, recursive = TRUE), add = TRUE)

  in_smi  <- file.path(workdir, "input.smi")
  out_csv <- file.path(workdir, "opera_out.csv")
  write_smi(compounds, in_smi)

  args <- c("-i", in_smi, "-o", out_csv,
            "-e", endpoints,
            "-v", as.character(verbose))
  if (isTRUE(neighbors)) args <- c(args, "-n")
  if (isTRUE(exp))       args <- c(args, "-x")

  res$command <- paste(shQuote(cmd), paste(shQuote(args), collapse = " "))

  log <- tryCatch(
    system2(cmd, args = args, stdout = TRUE, stderr = TRUE, timeout = timeout),
    error   = function(e) { res$error <<- conditionMessage(e); character(0) },
    warning = function(w) { res$error <<- conditionMessage(w); character(0) }
  )
  res$log <- as.character(log)

  if (!file.exists(out_csv)) {
    if (is.null(res$error)) {
      res$error <- "OPERA finished but produced no output file. See the run log."
    }
    return(res)
  }

  raw <- tryCatch(
    utils::read.csv(out_csv, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) { res$error <<- paste("Failed to read OPERA output:", conditionMessage(e)); NULL }
  )
  if (is.null(raw)) return(res)

  res$raw  <- raw
  res$tidy <- tryCatch(parse_opera_output(raw),
                       error = function(e) { res$error <<- conditionMessage(e); NULL })
  res$ok   <- !is.null(res$tidy)
  res
}

# ---------------------------------------------------------------------------
# Demo backend: clearly-labelled SIMULATED predictions so the UI can be
# exercised before OPERA is installed. These are NOT real OPERA results.
# ---------------------------------------------------------------------------
run_demo <- function(compounds, endpoints = ADME_ENDPOINTS) {
  set.seed(sum(utils::adist("opera")) + nrow(compounds))
  n <- nrow(compounds)
  ids <- if (!is.null(compounds$ID)) compounds$ID else paste0("Mol_", seq_len(n))

  mk <- function(tok, pred, lo, hi) {
    ad <- stats::runif(n)
    data.frame(
      row_id = seq_len(n), MoleculeID = ids, SMILES = compounds$SMILES,
      endpoint = tok,
      label = ADME_ENDPOINTS[[toupper(tok)]]$label %||% tok,
      unit  = ADME_ENDPOINTS[[toupper(tok)]]$unit  %||% "",
      pred = pred, predRange = sprintf("[%.2f, %.2f]", lo, hi),
      AD = as.numeric(ad > 0.25),
      AD_index = round(0.5 + 0.5 * ad, 2),
      Conf_index = round(0.4 + 0.5 * stats::runif(n), 2),
      stringsAsFactors = FALSE
    )
  }
  tidy <- rbind(
    mk("CACO2", round(stats::rnorm(n, -4.8, 0.4), 2),
       -5.6, -4.0),
    mk("FUB",   round(pmin(pmax(stats::rbeta(n, 2, 5), 0), 1), 3),
       0, 1),
    mk("Clint", round(abs(stats::rnorm(n, 25, 20)), 1),
       0, 120)
  )
  list(ok = TRUE, tidy = tidy, raw = tidy_to_wide(tidy),
       command = "DEMO MODE — simulated values (OPERA not called)",
       log = "Demo mode: predictions are randomly generated placeholders, NOT OPERA output.",
       error = NULL, demo = TRUE)
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a
