# OPERA ADME Predictor (Shiny)

A small R/Shiny app that runs the **ADME** models from
[NIEHS/OPERA](https://github.com/NIEHS/OPERA) on a set of structures (SMILES)
and presents the predictions with applicability-domain (AD) and confidence
indices.

## Does OPERA predict "ADME substrate and inhibition"?

**Not in the CYP450 sense.** OPERA has **no CYP450 substrate or inhibition
models** (no CYP1A2 / 2C9 / 2C19 / 2D6 / 3A4 flags). That kind of
substrate/inhibitor modelling comes from other resources such as
[NCATS ADME](https://opendata.ncats.nih.gov/adme), not from OPERA.

What OPERA **does** model for ADME are **three** endpoints — and that is
exactly what this app exposes:

| Endpoint | Meaning | Unit |
|----------|---------|------|
| **Caco-2** | Apparent permeability, logPapp (passive intestinal absorption proxy) | log(10⁻⁶ cm/s) |
| **FuB** | Fraction unbound in human plasma | fraction 0–1 |
| **Clint** | Human hepatic intrinsic clearance | µL/min/10⁶ cells |

(OPERA also predicts many physchem / environmental-fate / tox endpoints —
LogP, pKa, LogD, BP, MP, VP, WS, HL, KOA, Koc, BCF, biodegradation, CERAPP,
CoMPARA, CATMoS, etc. — but this app is deliberately scoped to the ADME three.)

## How it works

The app shells out to a locally installed **OPERA command-line executable**:

```
OPERA -i input.smi -o out.csv -e CACO2 FUB Clint -v 1
```

It then parses OPERA's output CSV. Output column names drift between OPERA
versions (e.g. `CACO2_pred`, `AD_CACO2`, `Conf_index_CACO2`), so the parser
**discovers** columns by pattern rather than hard-coding exact names — it
works across OPERA 2.6–2.9.

## Requirements

- **R ≥ 4.1** (uses the native `|>` pipe)
- R packages: `shiny`, `bslib`, `DT`, `bsicons`
- **OPERA** installed and on your `PATH` (or point the app at its executable)

Install the R packages:

```r
install.packages(c("shiny", "bslib", "DT", "bsicons"))
```

Install OPERA from the [official releases](https://github.com/NIEHS/OPERA) —
it ships a command-line installer for Windows and Linux (it bundles the
MATLAB Runtime).

## Running the app

```r
shiny::runApp(".")
```

or from a shell:

```bash
Rscript -e 'shiny::runApp(".", launch.browser = TRUE)'
```

### Pointing the app at OPERA

The app resolves the OPERA command in this order:

1. The **OPERA executable** field in the sidebar (under *OPERA backend*)
2. The `OPERA_CMD` environment variable
3. `OPERA` on your `PATH`

```bash
export OPERA_CMD=/opt/OPERA/application/OPERA   # example
```

### Demo mode

If OPERA is not installed yet, tick **Demo mode** in the sidebar to exercise
the interface with **simulated** values. These are random placeholders, clearly
labelled, and are **not** real OPERA predictions.

## Input formats

- **Paste SMILES**: one per line; an optional name may follow a comma or tab,
  e.g. `CC(=O)OC1=CC=CC=C1C(=O)O, Aspirin`
- **Upload**: `.smi` / `.txt` (SMILES per line) or `.csv` with a `SMILES` column
  (and optionally an `ID`/`Name` column). See [`examples/example.smi`](examples/example.smi).

## Project layout

```
opera_shiny/
├── app.R                # UI + server
├── R/
│   ├── adme_meta.R      # the 3 ADME endpoints (labels, units, interpretation)
│   ├── opera_cli.R      # runs the OPERA executable (+ demo backend)
│   ├── parse_opera.R    # version-tolerant output parser
│   └── helpers.R        # input parsing + UI helpers
├── examples/example.smi
└── README.md
```

## Caveats

- Predictions are only as good as OPERA's models; always check the
  **applicability domain** and **confidence index**.
- Caco-2 units/sign convention can vary by OPERA build — consult the OPERA
  QMRF for your installed version before interpreting absolute values.
