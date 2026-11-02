# =============================================================================
# OPERA ADME  -  Shiny front-end for the NIEHS/OPERA ADME models
#
# Predicts the THREE ADME endpoints that OPERA actually models:
#   * Caco-2 permeability (logPapp)
#   * Fraction unbound in plasma (FuB)
#   * Hepatic intrinsic clearance (Clint)
#
# NOTE: OPERA does NOT predict CYP450 substrate/inhibition. See README.
#
# Backend: shells out to a locally installed OPERA command-line executable.
# A clearly-labelled demo mode generates simulated values for UI review only.
# =============================================================================

library(shiny)
library(bslib)
library(DT)

# Source modules (works both when run from this dir and via runApp()).
for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f, local = FALSE)

EXAMPLE_SMILES <- paste(
  "CN1C=NC2=C1C(=O)N(C(=O)N2C)C\tCaffeine",
  "CC(=O)OC1=CC=CC=C1C(=O)O\tAspirin",
  "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O\tIbuprofen",
  "CC(=O)NC1=CC=C(C=C1)O\tParacetamol",
  "C1=CC=C(C=C1)C2=CC(=O)C3=C(C=C(C=C3O2)O)O\tChrysin",
  sep = "\n"
)

# ---- UI --------------------------------------------------------------------
ui <- page_sidebar(
  title = tagList(bsicons::bs_icon("capsule-pill"), "OPERA ADME Predictor"),
  # Bootstrap 5 default (system) fonts only -- no Google Fonts fetch, so the
  # app renders in offline / restricted-network environments.
  theme = bs_theme(version = 5, primary = "#2C6E8F"),
  window_title = "OPERA ADME Predictor",

  sidebar = sidebar(
    width = 360,
    title = "Input",

    radioButtons("input_mode", NULL,
                 choices = c("Paste SMILES" = "text", "Upload file" = "file"),
                 inline = TRUE),

    conditionalPanel(
      "input.input_mode == 'text'",
      textAreaInput("smiles_text", NULL, value = EXAMPLE_SMILES, rows = 8,
                    placeholder = "One SMILES per line, optional name after a comma/tab"),
      actionLink("load_example", "Load example set")
    ),
    conditionalPanel(
      "input.input_mode == 'file'",
      fileInput("smiles_file", NULL, accept = c(".smi", ".txt", ".csv"),
                buttonLabel = "Browse…"),
      helpText("Accepts .smi / .txt (SMILES per line) or .csv (a SMILES column).")
    ),

    hr(),
    accordion(
      open = FALSE,
      accordion_panel(
        "OPERA backend", icon = bsicons::bs_icon("gear"),
        textInput("opera_cmd", "OPERA executable",
                  value = Sys.getenv("OPERA_CMD", unset = "OPERA"),
                  placeholder = "e.g. /opt/OPERA/application/OPERA"),
        div(class = "d-flex gap-3",
            checkboxInput("opt_neighbors", "Nearest neighbours (-n)", FALSE),
            checkboxInput("opt_exp", "Experimental (-x)", FALSE)),
        sliderInput("opt_verbose", "Verbosity (-v)", 0, 2, 1, step = 1),
        checkboxInput("demo_mode", "Demo mode (simulated — no OPERA)", FALSE),
        uiOutput("backend_status")
      )
    ),

    hr(),
    input_task_button("run", "Predict ADME", icon = bsicons::bs_icon("play-fill")),
    uiOutput("n_compounds")
  ),

  # Factual banner about OPERA's ADME scope.
  div(class = "alert alert-info d-flex align-items-center py-2 mb-3",
      bsicons::bs_icon("info-circle", class = "me-2"),
      HTML(paste0(
        "OPERA models <b>three</b> ADME endpoints — Caco-2 permeability, ",
        "fraction unbound (FuB) and intrinsic clearance (Clint). ",
        "It does <b>not</b> predict CYP450 substrate or inhibition."))),

  navset_card_tab(
    id = "tabs",
    nav_panel(
      "Overview", icon = bsicons::bs_icon("grid-1x2"),
      uiOutput("summary_boxes"),
      card(
        card_header("Predictions"),
        DTOutput("results_tbl"),
        card_footer(downloadButton("dl_csv", "Download results (CSV)",
                                    class = "btn-sm btn-primary"))
      )
    ),
    nav_panel(
      "Endpoint detail", icon = bsicons::bs_icon("clipboard-data"),
      uiOutput("compound_picker"),
      uiOutput("detail_cards")
    ),
    nav_panel(
      "Run log", icon = bsicons::bs_icon("terminal"),
      card(card_header("Command"), verbatimTextOutput("cmd_out")),
      card(card_header("OPERA output / messages"), verbatimTextOutput("log_out"))
    ),
    nav_panel(
      "About", icon = bsicons::bs_icon("book"),
      card(card_body(uiOutput("about_md")))
    )
  )
)

# ---- Server ----------------------------------------------------------------
server <- function(input, output, session) {

  observeEvent(input$load_example, {
    updateTextAreaInput(session, "smiles_text", value = EXAMPLE_SMILES)
  })

  compounds <- reactive({
    if (identical(input$input_mode, "file")) {
      req(input$smiles_file)
      parse_smiles_file(input$smiles_file$datapath, input$smiles_file$name)
    } else {
      parse_smiles_text(input$smiles_text)
    }
  })

  output$n_compounds <- renderUI({
    n <- nrow(compounds())
    helpText(sprintf("%d valid compound%s detected.", n, if (n == 1) "" else "s"))
  })

  output$backend_status <- renderUI({
    if (isTRUE(input$demo_mode))
      return(span(class = "badge bg-warning text-dark", "Demo mode active"))
    if (opera_available(input$opera_cmd))
      span(class = "badge bg-success", "OPERA executable found")
    else
      span(class = "badge bg-danger", "OPERA not found on this host")
  })

  # Run predictions on button press.
  results <- eventReactive(input$run, {
    cmp <- compounds()
    validate(need(nrow(cmp) > 0, "Enter at least one SMILES string."))
    if (isTRUE(input$demo_mode)) {
      run_demo(cmp)
    } else {
      run_opera(cmp,
                cmd       = input$opera_cmd,
                neighbors = input$opt_neighbors,
                exp       = input$opt_exp,
                verbose   = as.integer(input$opt_verbose))
    }
  })

  observeEvent(results(), {
    r <- results()
    if (!isTRUE(r$ok)) {
      showNotification(paste("Prediction failed:", r$error %||% "unknown error"),
                       type = "error", duration = 8)
    } else if (isTRUE(r$demo)) {
      showNotification("Demo mode: values are simulated, not OPERA predictions.",
                       type = "warning", duration = 6)
    } else {
      showNotification("OPERA predictions complete.", type = "message", duration = 4)
    }
  })

  tidy <- reactive({ r <- results(); req(isTRUE(r$ok)); r$tidy })

  # --- summary value boxes (mean per endpoint) ---
  output$summary_boxes <- renderUI({
    t <- tidy()
    boxes <- lapply(ADME_ENDPOINTS, function(ep) {
      sub <- t[t$endpoint == ep$token, ]
      if (!nrow(sub)) return(NULL)
      n_in <- sum(sub$AD >= 0.5, na.rm = TRUE)
      val  <- if (nrow(sub) == 1) fmt_num(sub$pred[1]) else fmt_num(mean(sub$pred, na.rm = TRUE))
      value_box(
        title = ep$label,
        value = val,
        showcase = bsicons::bs_icon(ep$icon),
        theme = "primary",
        p(class = "small mb-0", ep$unit),
        p(class = "small mb-0", sprintf("%d / %d in applicability domain", n_in, nrow(sub)))
      )
    })
    do.call(layout_column_wrap,
            c(list(width = 1/3), Filter(Negate(is.null), boxes)))
  })

  # --- results table (wide) ---
  output$results_tbl <- renderDT({
    w <- tidy_to_wide(tidy())
    datatable(w, rownames = FALSE, selection = "single",
              options = list(pageLength = 10, scrollX = TRUE, dom = "tip")) |>
      formatStyle(columns = grep("_AD$", names(w), value = TRUE),
                  backgroundColor = styleEqual(c("in", "out"),
                                               c("#d1e7dd", "#fff3cd")))
  })

  # --- compound picker for detail tab ---
  output$compound_picker <- renderUI({
    t <- tidy()
    ids <- unique(t$MoleculeID)
    selectInput("detail_id", "Compound", choices = ids, width = "360px")
  })

  output$detail_cards <- renderUI({
    t <- tidy(); req(input$detail_id)
    sub <- t[t$MoleculeID == input$detail_id, ]
    cards <- lapply(split(sub, sub$endpoint), function(row) {
      ep <- ADME_ENDPOINTS[[toupper(row$endpoint[1])]]
      if (is.null(ep)) ep <- list(label = row$label[1], unit = row$unit[1],
                                  detail = "", interpret = function(v) NA_character_,
                                  icon = "clipboard-data")
      note <- tryCatch(ep$interpret(row$pred[1]), error = function(e) NA_character_)
      card(
        card_header(tagList(bsicons::bs_icon(ep$icon), " ", ep$label)),
        card_body(
          div(class = "d-flex justify-content-between align-items-baseline",
              div(h2(class = "mb-0", fmt_num(row$pred[1])),
                  span(class = "text-muted small", ep$unit)),
              ad_badge(row$AD[1])),
          if (!is.na(row$predRange[1]) && nzchar(row$predRange[1]))
            p(class = "small text-muted mb-1", "Predicted range: ", row$predRange[1]),
          if (!is.na(note)) div(class = "alert alert-light py-1 px-2 small mb-2", note),
          div(class = "small mb-1",
              "Confidence index: ",
              span(class = paste("badge", conf_class(row$Conf_index[1])),
                   fmt_num(row$Conf_index[1], 2)),
              if (!is.na(row$AD_index[1]))
                tagList(" · AD index: ", span(class = "badge bg-secondary",
                                                   fmt_num(row$AD_index[1], 2)))),
          p(class = "text-muted small mb-0", ep$detail)
        )
      )
    })
    do.call(layout_column_wrap, c(list(width = 1/3), cards))
  })

  # --- run log ---
  output$cmd_out <- renderText({ r <- results(); r$command %||% "" })
  output$log_out <- renderText({
    r <- results()
    paste(c(if (!is.null(r$error)) paste("ERROR:", r$error), r$log), collapse = "\n")
  })

  # --- download ---
  output$dl_csv <- downloadHandler(
    filename = function() sprintf("opera_adme_%s.csv", Sys.Date()),
    content = function(file) utils::write.csv(tidy_to_wide(tidy()), file, row.names = FALSE)
  )

  # --- about ---
  output$about_md <- renderUI({
    HTML(paste0(
      "<h4>About this app</h4>",
      "<p>This app is a front-end for <a href='https://github.com/NIEHS/OPERA' target='_blank'>",
      "NIEHS/OPERA</a>, running its three ADME QSAR models and presenting the ",
      "predictions with applicability-domain (AD) and confidence indices.</p>",
      "<h5>The three OPERA ADME endpoints</h5><ul>",
      paste(vapply(ADME_ENDPOINTS, function(ep)
        sprintf("<li><b>%s</b> (%s) — %s</li>", ep$label, ep$unit, ep$desc),
        character(1)), collapse = ""),
      "</ul>",
      "<h5>What OPERA does <i>not</i> do</h5>",
      "<p>OPERA has <b>no CYP450 substrate or inhibition models</b> ",
      "(CYP1A2/2C9/2C19/2D6/3A4). Such models exist elsewhere — e.g. ",
      "<a href='https://opendata.ncats.nih.gov/adme' target='_blank'>NCATS ADME</a> ",
      "— but they are not part of OPERA.</p>",
      "<h5>Interpreting AD &amp; confidence</h5>",
      "<p><b>In domain</b> means the compound resembles OPERA's training set; ",
      "predictions outside the domain are less reliable. The <b>confidence index</b> ",
      "(0–1, higher = better) combines AD and local model performance.</p>"
    ))
  })
}

shinyApp(ui, server)
