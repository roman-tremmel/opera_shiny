# parse_opera.R -------------------------------------------------------------
# Turn a raw OPERA output CSV (data.frame) into a tidy long table, one row per
# (compound, endpoint). OPERA column names drift between versions
# (e.g. "CACO2_pred", "AD_CACO2", "AD_index_CACO2", "Conf_index_CACO2"), so we
# DISCOVER columns by pattern instead of hard-coding exact strings.
# ---------------------------------------------------------------------------

# Find the single best-matching column for one "field" of an endpoint.
# `cols`  : available column names
# `token` : endpoint token, e.g. "CACO2" (matched case-insensitively)
# `field` : one of "pred", "range", "ad", "ad_index", "conf"
.find_col <- function(cols, token, field) {
  tok <- tolower(token)
  lc  <- tolower(cols)
  has <- grepl(tok, lc, fixed = TRUE)
  if (!any(has)) return(NA_character_)
  cand <- cols[has]
  lcc  <- lc[has]

  pick <- switch(
    field,
    "conf"     = grepl("conf", lcc),
    "ad_index" = grepl("ad[_.]?index", lcc),
    "ad"       = grepl("^ad[_.]", lcc) & !grepl("index", lcc),
    "range"    = grepl("range", lcc),
    "pred"     = grepl("pred", lcc) & !grepl("range", lcc) &
                 !grepl("ad[_.]?index|conf", lcc),
    logical(length(lcc))
  )

  if (any(pick)) return(cand[which(pick)[1]])

  # Fallback for the predicted value: a bare column equal to the token
  # (some builds emit just "FUB" rather than "FUB_pred").
  if (field == "pred") {
    bare <- lcc == tok
    if (any(bare)) return(cand[which(bare)[1]])
  }
  NA_character_
}

# Identify the molecule-id / SMILES columns present in OPERA output.
.id_cols <- function(cols) {
  lc <- tolower(cols)
  list(
    id    = cols[which(grepl("molecule.?id|^id$|^name$", lc))[1]],
    smi   = cols[which(grepl("original.*smiles|input.*smiles|^smiles$", lc))[1]],
    canon = cols[which(grepl("qsar.?ready.*smiles|canonical.*smiles", lc))[1]]
  )
}

#' Parse an OPERA output data.frame into a tidy long table.
#'
#' @param df         data.frame read from OPERA's output CSV
#' @param endpoints  named list like ADME_ENDPOINTS
#' @return data.frame with columns:
#'   row_id, MoleculeID, SMILES, endpoint (token), label, unit,
#'   pred, predRange, AD, AD_index, Conf_index
parse_opera_output <- function(df, endpoints = ADME_ENDPOINTS) {
  stopifnot(is.data.frame(df))
  cols <- names(df)
  ids  <- .id_cols(cols)

  mol_id <- if (!is.na(ids$id)) as.character(df[[ids$id]]) else as.character(seq_len(nrow(df)))
  smiles <- if (!is.na(ids$smi)) as.character(df[[ids$smi]])
            else if (!is.na(ids$canon)) as.character(df[[ids$canon]])
            else NA_character_

  out <- list()
  for (ep in endpoints) {
    c_pred  <- .find_col(cols, ep$token, "pred")
    c_range <- .find_col(cols, ep$token, "range")
    c_ad    <- .find_col(cols, ep$token, "ad")
    c_adi   <- .find_col(cols, ep$token, "ad_index")
    c_conf  <- .find_col(cols, ep$token, "conf")

    # Skip endpoints OPERA did not emit at all.
    if (is.na(c_pred) && is.na(c_ad) && is.na(c_conf)) next

    num <- function(cn) if (is.na(cn)) rep(NA_real_, nrow(df)) else suppressWarnings(as.numeric(df[[cn]]))

    out[[ep$token]] <- data.frame(
      row_id     = seq_len(nrow(df)),
      MoleculeID = mol_id,
      SMILES     = smiles,
      endpoint   = ep$token,
      label      = ep$label,
      unit       = ep$unit,
      pred       = num(c_pred),
      predRange  = if (is.na(c_range)) NA_character_ else as.character(df[[c_range]]),
      AD         = num(c_ad),
      AD_index   = num(c_adi),
      Conf_index = num(c_conf),
      stringsAsFactors = FALSE
    )
  }

  if (!length(out)) {
    stop("No ADME endpoint columns (Caco-2 / FuB / Clint) were found in the OPERA output. ",
         "Columns present: ", paste(cols, collapse = ", "))
  }
  do.call(rbind, out)
}

# Wide (one row per compound) view built from the tidy table, for a results
# grid: MoleculeID + <token>_pred / <token>_AD columns.
tidy_to_wide <- function(tidy) {
  base <- unique(tidy[, c("row_id", "MoleculeID", "SMILES")])
  base <- base[order(base$row_id), , drop = FALSE]
  for (tok in unique(tidy$endpoint)) {
    sub <- tidy[tidy$endpoint == tok, ]
    m   <- match(base$row_id, sub$row_id)
    base[[paste0(tok, "_pred")]] <- round(sub$pred[m], 3)
    base[[paste0(tok, "_AD")]]   <- ifelse(is.na(sub$AD[m]), NA,
                                           ifelse(sub$AD[m] >= 0.5, "in", "out"))
    base[[paste0(tok, "_Conf")]] <- round(sub$Conf_index[m], 2)
  }
  base$row_id <- NULL
  base
}
