# adme_meta.R ---------------------------------------------------------------
# Metadata for the three ADME endpoints that NIEHS/OPERA actually predicts.
#
# OPERA does NOT model CYP450 substrate/inhibition. Its released ADME models
# are exactly these three:
#   * Caco-2 apparent permeability (logPapp)
#   * Fraction unbound in plasma  (FuB)
#   * Hepatic intrinsic clearance (Clint)
# ---------------------------------------------------------------------------

# Each entry:
#   token    - case-insensitive string used both as the OPERA `-e` endpoint
#              name and to discover the endpoint's columns in OPERA output.
#   label    - human-readable name
#   short    - compact label for tables / value boxes
#   unit     - measurement unit as reported by OPERA
#   icon     - bsicons name
#   desc     - one-line description
#   detail   - longer explanation shown in the "about" panel
#   interpret- function(value) -> short qualitative note (character) or NA
ADME_ENDPOINTS <- list(
  CACO2 = list(
    token  = "CACO2",
    label  = "Caco-2 Permeability",
    short  = "Caco-2",
    unit   = "log(10⁻⁶ cm/s)",
    icon   = "arrow-left-right",
    desc   = "Apparent permeability (logPapp) across a Caco-2 cell monolayer — a proxy for passive intestinal absorption.",
    detail = paste(
      "OPERA predicts log Papp. Higher values indicate greater predicted",
      "passive permeability and, generally, better intestinal absorption.",
      "Exact units and sign convention depend on the OPERA model version —",
      "consult the OPERA QMRF for the installed build."
    ),
    interpret = function(v) {
      if (is.na(v)) return(NA_character_)
      # Soft, directional guidance only (thresholds vary by version).
      if (v >= -4.80) "Relatively high predicted permeability"
      else if (v >= -5.15) "Moderate predicted permeability"
      else "Relatively low predicted permeability"
    }
  ),

  FUB = list(
    token  = "FUB",
    label  = "Fraction Unbound in Plasma",
    short  = "FuB",
    unit   = "fraction (0–1)",
    icon   = "droplet-half",
    desc   = "Fraction of compound circulating free (unbound to plasma proteins) in human plasma.",
    detail = paste(
      "FuB ranges from 0 to 1. Values near 0 indicate extensive plasma-protein",
      "binding (little free drug); values near 1 indicate the compound is",
      "largely unbound. FuB is a key input for IVIVE and PBPK workflows."
    ),
    interpret = function(v) {
      if (is.na(v)) return(NA_character_)
      if (v < 0.05) "Highly protein-bound (very little free fraction)"
      else if (v < 0.30) "Substantially protein-bound"
      else if (v < 0.70) "Moderately bound"
      else "Largely unbound / free"
    }
  ),

  CLINT = list(
    token  = "Clint",
    label  = "Hepatic Intrinsic Clearance",
    short  = "Clint",
    unit   = "µL/min/10⁶ cells",
    icon   = "funnel",
    desc   = "Human hepatic intrinsic clearance (hepatocyte) — the capacity of the liver to metabolise the compound.",
    detail = paste(
      "Clint reflects metabolic turnover by hepatocytes. Higher values imply",
      "faster hepatic metabolism and typically lower metabolic stability.",
      "Reported in µL/min per 10⁶ cells."
    ),
    interpret = function(v) {
      if (is.na(v)) return(NA_character_)
      if (v < 5) "Low clearance / metabolically stable"
      else if (v < 50) "Moderate clearance"
      else "High clearance / rapidly metabolised"
    }
  )
)

# Ordered token vector passed to OPERA's `-e` flag.
adme_tokens <- function() vapply(ADME_ENDPOINTS, function(x) x$token, character(1))
