library(shiny)
library(shinyWidgets)
library(data.table)
library(dplyr)
library(tidyr)          # fix: explicit import (used for pivot_longer)
library(plotly)
library(DT)

## ---------------------------------------------------------------
## Data paths per ancestry
## ---------------------------------------------------------------
ANCESTRY_FILES <- list(
  "European"  = list(ld = "../results/white/UKB_white_ld_residuals.txt",
                     gw = "../results/white/lpa_man_white.regenie.gz"),
  "African"   = list(ld = "../results/african/UKB_afr_ld_residuals.txt",
                     gw = "../results/african/lpa_man.regenie_afr.gz"),
  "Asian"     = list(ld = "../results/asian/UKB_asian_ld_residuals.txt",
                     gw = "../results/asian/lpa_man.regenie_asian.gz")
)

## ---------------------------------------------------------------
## Colours (constant across ancestries)
## ---------------------------------------------------------------
COL_SELECTED <- "#f39c12"   # orange — selected from table

# Ancestry colours: VNTR (light) and non-rep (dark) tones
ANCESTRY_COLS <- list(
  European = list(vntr = "#aed6f1", nonrep = "#2980b9"),
  African  = list(vntr = "#f0b27a", nonrep = "#e67e22"),
  Asian    = list(vntr = "#a9dfbf", nonrep = "#27ae60")
)

## fix: extract duplicated VNTR position logic into one place
is_vntr_pos <- function(pos) {
  !is.na(pos) & ((pos >= 480 & pos <= 840) | (pos >= 4643 & pos <= 5025))
}

## Helper: load and process data for a given ancestry
load_ancestry <- function(ancestry) {
  files     <- ANCESTRY_FILES[[ancestry]]
  ld_mat    <- as.matrix(fread(files$ld, header = TRUE))
  regenie   <- fread(files$gw, header = TRUE)

  # Use RSID where ID is missing (".")
  if ("RSID" %in% names(regenie))
    regenie[ID == ".", ID := RSID]

  rownames(ld_mat) <- colnames(ld_mat) <- regenie$ID
  ld_r2     <- ld_mat^2
  all_ids   <- rownames(ld_r2)
  id_pos    <- suppressWarnings(as.integer(sub("^6:([0-9]+)[_:].*", "\\1", all_ids)))
  is_vntr   <- is_vntr_pos(id_pos)
  snp_types <- ifelse(is_vntr, "VNTR region", "Non-repetitive region")
  vntr_ids  <- all_ids[is_vntr][order(id_pos[is_vntr])]

  snp_meta  <- regenie %>%
    select(ID, GENPOS, ALLELE0, ALLELE1, A1FREQ, BETA, SE, LOG10P) %>%
    mutate(
      id_pos = suppressWarnings(as.integer(sub("^6:([0-9]+)[_:].*", "\\1", ID))),
      Type   = ifelse(is_vntr_pos(id_pos), "VNTR region", "Non-repetitive region")
    ) %>%
    select(-id_pos) %>%
    distinct(ID, .keep_all = TRUE)

  snp_type_map <- setNames(snp_meta$Type, snp_meta$ID)

  list(ld_r2 = ld_r2, all_ids = all_ids, is_vntr = is_vntr,
       snp_types = snp_types, vntr_ids = vntr_ids,
       snp_meta = snp_meta, snp_type_map = snp_type_map)
}

## fix: pre-load all ancestries once at startup — no redundant re-reads
ANCESTRY_DATA <- lapply(names(ANCESTRY_FILES), load_ancestry)
names(ANCESTRY_DATA) <- names(ANCESTRY_FILES)

## ---------------------------------------------------------------
## UI
## ---------------------------------------------------------------
ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-family: 'Helvetica Neue', sans-serif; }
    .well { background-color: #f8f9fa; border: none; }
    h4 { color: #2c3e50; }
    .legend-box {
      display: inline-block; width: 12px; height: 12px;
      border-radius: 2px; margin-right: 5px;
    }
  "))),

  tags$div(
    style = "padding: 12px 0 8px 0;",
    tags$img(src = "logo.svg", height = "55px", alt = "Genetic Architecture of LPA")
  ),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("Query SNP"),
      p("Query a SNP to explore its LD partners."),
      selectizeInput("query_snp", "VNTR SNP",
        choices  = ANCESTRY_DATA[["European"]]$vntr_ids,  # fix: use pre-loaded data
        selected = NULL,
        options  = list(
          placeholder    = "Search / select VNTR SNP …",
          onInitialize   = I('function() { this.setValue(""); }')
        )
      ),
      tags$p(style = "text-align:center; color:#888; margin: 2px 0;", "— or —"),
      textInput("rs_input", "rs number", placeholder = "e.g. rs78822335"),

      hr(),

      h4("Select Ancestry"),
      radioGroupButtons("ancestry", NULL,
        choices   = names(ANCESTRY_FILES),
        selected  = "European",
        direction = "horizontal",
        size      = "sm",
        checkIcon = list(yes = icon("check")),
        choiceNames = list(
          tags$span(style = paste0("color:", ANCESTRY_COLS$European$nonrep, "; font-weight:600;"), "European"),
          tags$span(style = paste0("color:", ANCESTRY_COLS$African$nonrep,  "; font-weight:600;"), "African"),
          tags$span(style = paste0("color:", ANCESTRY_COLS$Asian$nonrep,    "; font-weight:600;"), "Asian")
        ),
        choiceValues = names(ANCESTRY_FILES)
      ),

      conditionalPanel(
        condition = "input.tabs == 'LD Heatmap'",
        hr(),
        h4("Display filters"),
        checkboxGroupInput("show_types", "SNP types in heatmap",
          choices  = c("VNTR region", "Non-repetitive region"),
          selected = c("VNTR region", "Non-repetitive region")
        )
      ),

      hr(),

      h4("Selected SNP info"),
      tableOutput("query_info")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        ## --- Tab 1: LD partners ---
        tabPanel("LD Partners",
          br(),
          fluidRow(
            column(4,
              sliderInput("min_r2", "Min R² to display",
                min = 0, max = 1, value = 0.1, step = 0.05)
            ),
            column(4,
              checkboxGroupInput("partner_types", "Show partner types",
                choices  = c("VNTR region", "Non-repetitive region"),
                selected = c("VNTR region", "Non-repetitive region")
              )
            )
            # fix: removed static HTML legend — bar plot plotly legend is correct
            #      and updates with ancestry; the old static one used hardcoded
            #      European colours
          ),
          plotlyOutput("bar_plot", height = "350px"),
          br(),
          downloadButton("download_csv", "Download CSV", class = "btn-sm"),
          br(), br(),
          DTOutput("ld_table")
        ),

        ## --- Tab 2: Heatmap ---
        tabPanel("LD Heatmap",
          br(),
          fluidRow(
            column(4,
              sliderInput("r2_highlight", "Highlight rs SNPs in LD ≥ R²",
                min = 0, max = 1, value = 0.7, step = 0.05)
            )
          ),
          plotlyOutput("heatmap", height = "680px")
        ),

        ## --- Tab 3: About ---
        tabPanel("About",
          br(),
          tags$div(style = "max-width: 760px;",

            h4(tags$a("Institute of Genetic Epidemiology, Medical University of Innsbruck",
                href = "https://genepi.i-med.ac.at", target = "_blank")),

            tags$p(tags$strong("Silvia Di Maio,"),
              tags$a("silvia.di-maio@i-med.ac.at", href = "mailto:silvia.di-maio@i-med.ac.at")),
            tags$p(tags$strong("Johanna F. Schachtl-Rieß,"),
              tags$a("johanna.schachtl-riess@i-med.ac.at", href = "mailto:johanna.schachtl-riess@i-med.ac.at")),
            tags$p(tags$strong("Stefan Coassin,"),
              tags$a("stefan.coassin@i-med.ac.at", href = "mailto:stefan.coassin@i-med.ac.at")),
            tags$p(tags$strong("Sebastian Schönherr,"),
              tags$a("sebastian.schoenherr@i-med.ac.at", href = "mailto:sebastian.schoenherr@i-med.ac.at")),

            hr(),

            h4("Data"),
            tags$p(
              "GWAS summary statistics were derived from",
              tags$a("UK Biobank", href = "https://www.ukbiobank.ac.uk", target = "_blank"),
              "data (application number 62905). Analysis code and pipeline are publicly available on",
              tags$a("GitHub.", href = "https://github.com/seppinho/vntr-gwas-framework/", target = "_blank")
            ),

            hr(),

            h4("Funding"),
            tags$p(
              "This work is supported by the Austrian Science Fund (FWF) under grant",
              tags$a("PAT3357425", href = "https://www.fwf.ac.at/forschungsradar/10.55776/PAT3357425", target = "_blank"),
              "—", em("Intra-Repeat VNTR Resolution and Risk Prediction"),
              "(PI: Sebastian Schönherr, Medical University of Innsbruck)."
            )
          )
        )
      )
    )
  )
)

## ---------------------------------------------------------------
## Server
## ---------------------------------------------------------------
server <- function(input, output, session) {

  ## fix: reactiveVals moved to top of server for clarity
  query_snp  <- reactiveVal(NULL)
  sel_snp_id <- reactiveVal(NULL)   # ID selected from the LD table

  ## -- Reactive: data for selected ancestry (now just a lookup) --
  adata <- reactive({
    ANCESTRY_DATA[[input$ancestry]]
  })

  # Update inputs when ancestry changes, preserving query SNP if it exists
  observeEvent(input$ancestry, {
    d       <- adata()
    current <- query_snp()
    if (!is.null(current) && current %in% d$all_ids) {
      # SNP exists in new ancestry — keep it
      updateSelectizeInput(session, "query_snp", choices = d$vntr_ids,
                           selected = if (current %in% d$vntr_ids) current else character(0),
                           options = list(placeholder  = "Search / select VNTR SNP …",
                                          onInitialize = I('function() { this.setValue(""); }')))
    } else {
      # SNP not in new ancestry — reset
      query_snp(NULL)
      updateSelectizeInput(session, "query_snp", choices = d$vntr_ids,
                           selected = character(0),
                           options = list(placeholder  = "Search / select VNTR SNP …",
                                          onInitialize = I('function() { this.setValue(""); }')))
      updateTextInput(session, "rs_input", value = "")
    }
  })

  ## -- Reactive: which SNPs to show in the heatmap --
  vis_ids <- reactive({
    d <- adata()
    d$all_ids[d$snp_types %in% input$show_types]
  })

  ## -- Reactive: sub-matrix for heatmap --
  sub_r2 <- reactive({
    s <- vis_ids()
    req(length(s) > 1)   # character(0) subscript crashes matrix indexing in R
    adata()$ld_r2[s, s]
  })

  ## -- Build heatmap --
  output$heatmap <- renderPlotly({
    mat  <- sub_r2()
    ids  <- rownames(mat)
    req(length(ids) > 1)

    types_vec <- adata()$snp_type_map[ids]

    # Per-cell hover text matrix: "row_id (type) | col_id (type)"
    row_label <- paste0(ids, "<br>", types_vec)
    col_label <- paste0(ids, "<br>", types_vec)
    hover_mat <- outer(row_label, col_label,
      FUN = function(r, c) paste0(r, "<br>vs<br>", c))

    # fix: removed unused tick_labels computation

    # Build highlight shapes based on selected query SNP
    shapes <- list()
    s <- query_snp()

    if (!is.null(s) && s %in% ids) {
      q_idx <- which(ids == s) - 1L   # 0-based index for plotly shapes

      # Outline query SNP column (red)
      shapes <- c(shapes, list(list(
        type      = "rect", xref = "x", yref = "paper",
        x0        = q_idx - 0.5, x1 = q_idx + 0.5,
        y0        = 0, y1 = 1,
        fillcolor = "rgba(231,76,60,0.08)",
        line      = list(color = "rgba(231,76,60,0.9)", width = 2)
      )))

      # Outline query SNP row (red)
      shapes <- c(shapes, list(list(
        type      = "rect", xref = "paper", yref = "y",
        x0        = 0, x1 = 1,
        y0        = q_idx - 0.5, y1 = q_idx + 0.5,
        fillcolor = "rgba(231,76,60,0.08)",
        line      = list(color = "rgba(231,76,60,0.9)", width = 2)
      )))

      # Blue bands for rs SNPs in LD with query >= threshold
      ld_row      <- adata()$ld_r2[s, ids]
      partner_ids <- names(ld_row)[
        ld_row >= input$r2_highlight &
        names(ld_row) != s &
        adata()$snp_type_map[names(ld_row)] == "Non-repetitive region"
      ]

      for (pid in partner_ids) {
        p_idx <- which(ids == pid) - 1L

        # Column band only (no row bands — keeps the heatmap readable)
        shapes <- c(shapes, list(list(
          type      = "rect", xref = "x", yref = "paper",
          x0        = p_idx - 0.5, x1 = p_idx + 0.5,
          y0        = 0, y1 = 1,
          fillcolor = "rgba(41,128,185,0.22)",
          line      = list(color = "rgba(41,128,185,0.9)", width = 1.5)
        )))
      }
    }

    plot_ly(
      x         = ids,
      y         = ids,
      z         = mat,
      text      = hover_mat,
      type      = "heatmap",
      colorscale = list(
        c(0,   "white"),
        c(0.3, "#ffffb2"),
        c(0.6, "#fd8d3c"),
        c(1,   "#bd0026")
      ),
      zmin = input$r2_highlight, zmax = 1,
      hovertemplate = paste0(
        "%{text}<br>",
        "<b>R² = %{z:.3f}</b>",
        "<extra></extra>"
      ),
      colorbar = list(title = "R²", len = 0.6)
    ) %>%
      layout(
        xaxis = list(
          showticklabels = FALSE,
          title          = "",
          showgrid       = FALSE
        ),
        yaxis = list(
          showticklabels = FALSE,
          title          = "",
          showgrid       = FALSE,
          autorange      = "reversed"
        ),
        shapes = shapes,
        margin = list(l = 160, b = 160, t = 30, r = 20)
      ) %>%
      event_register("plotly_click")
  })

  # Update from heatmap click
  observeEvent(event_data("plotly_click", source = "heatmap"), {
    click <- event_data("plotly_click", source = "heatmap")
    clicked_id <- click$x
    d <- adata()
    if (!is.null(clicked_id) && clicked_id %in% d$all_ids) {
      query_snp(clicked_id)
      if (clicked_id %in% d$vntr_ids) {   # fix: use d$vntr_ids directly
        updateSelectizeInput(session, "query_snp", selected = clicked_id)
        updateTextInput(session, "rs_input", value = "")
      } else {
        updateSelectizeInput(session, "query_snp", selected = character(0))
        updateTextInput(session, "rs_input", value = clicked_id)
      }
      updateTabsetPanel(session, "tabs", selected = "LD Partners")
    }
  })

  # Update from VNTR dropdown
  observeEvent(input$query_snp, ignoreInit = TRUE, {
    if (!is.null(input$query_snp) && input$query_snp != "") {
      query_snp(input$query_snp)
      updateTextInput(session, "rs_input", value = "")
    }
  })

  # Update from rs text input
  observeEvent(input$rs_input, ignoreInit = TRUE, {
    rs <- trimws(input$rs_input)
    if (nchar(rs) == 0) return()
    if (rs %in% adata()$all_ids) {
      query_snp(rs)
      updateSelectizeInput(session, "query_snp", selected = character(0))
      updateTabsetPanel(session, "tabs", selected = "LD Partners")
    } else {
      # fix: notify the user when the rs ID is not found
      showNotification(paste0('"', rs, '" not found in the current ancestry dataset.'),
                       type = "warning", duration = 4)
    }
  })

  ## -- Query SNP metadata --
  output$query_info <- renderTable({
    s <- query_snp()
    if (is.null(s)) return(data.frame(Info = "No SNP selected"))
    adata()$snp_meta %>%
      filter(ID == s) %>%
      select(ID, Type, GENPOS, ALLELE0, ALLELE1, A1FREQ, BETA, SE, LOG10P) %>%
      mutate(across(everything(), as.character)) %>%
      pivot_longer(everything(), names_to = "Field", values_to = "Value")
  }, striped = TRUE, hover = TRUE, bordered = TRUE, spacing = "xs")

  ## -- Reactive: all LD partners of query SNP (no UI filters) --
  partners_base <- reactive({
    s <- query_snp()
    d <- adata()
    req(s, s %in% rownames(d$ld_r2))

    row_vals <- d$ld_r2[s, ]
    row_vals <- row_vals[names(row_vals) != s]   # remove self

    data.frame(
      ID   = names(row_vals),
      r2   = as.numeric(row_vals),
      Type = d$snp_type_map[names(row_vals)]
    ) %>%
      left_join(d$snp_meta %>% select(ID, GENPOS, ALLELE0, ALLELE1,
                                      A1FREQ, BETA, SE, LOG10P),
                by = "ID")
  })

  ## -- Reactive: filtered view (reruns only when slider/checkboxes change) --
  partners_df <- reactive({
    partners_base() %>%
      filter(
        Type %in% input$partner_types,
        r2   >= input$min_r2
      ) %>%
      arrange(desc(r2))
  })

  ## fix: track selected row by ID rather than index — robust to DT filtering
  observeEvent(input$ld_table_rows_selected, {
    df      <- partners_df()
    sel_idx <- input$ld_table_rows_selected
    sel_snp_id(if (length(sel_idx) > 0) df[sel_idx, "ID"] else NULL)
  })

  ## -- Bar plot: top LD partners --
  output$bar_plot <- renderPlotly({
    df <- partners_df()
    req(nrow(df) > 0)

    df_top  <- head(df, 60)

    acols     <- ANCESTRY_COLS[[input$ancestry]]
    color_vec <- ifelse(df_top$Type == "VNTR region", acols$vntr, acols$nonrep)
    sel_id    <- sel_snp_id()
    if (!is.null(sel_id) && sel_id %in% df_top$ID)
      color_vec[df_top$ID == sel_id] <- COL_SELECTED

    plot_ly(
      data          = df_top,
      x             = ~reorder(ID, -r2),
      y             = ~r2,
      type          = "bar",
      marker        = list(
        color = color_vec,
        line  = list(color = "rgba(0,0,0,0.15)", width = 0.5)
      ),
      hovertemplate = ~paste0(
        "<b>", ID, "</b><br>",
        "Type: ", Type, "<br>",
        "Position: ", GENPOS, "<br>",
        "R² = ", round(r2, 3), "<br>",
        "log10P = ", round(LOG10P, 2),
        "<extra></extra>"
      ),
      showlegend = FALSE
    ) %>%
      # Dummy traces for legend — these use reactive ancestry colours, so the
      # legend stays in sync when the user switches ancestry (fix for old static
      # HTML legend that was hardcoded to European colours)
      add_trace(
        x = NA_character_, y = NA_real_, type = "bar",
        marker = list(color = acols$vntr),
        name = "VNTR region", showlegend = TRUE, inherit = FALSE
      ) %>%
      add_trace(
        x = NA_character_, y = NA_real_, type = "bar",
        marker = list(color = acols$nonrep),
        name = "Non-repetitive region", showlegend = TRUE, inherit = FALSE
      ) %>%
      add_trace(
        x = NA_character_, y = NA_real_, type = "bar",
        marker = list(color = COL_SELECTED),
        name = "Selected", showlegend = TRUE, inherit = FALSE
      ) %>%
      layout(
        title  = list(
          text = query_snp(),
          font = list(size = 14)
        ),
        xaxis   = list(title = "", tickangle = -50, tickfont = list(size = 8)),
        yaxis   = list(title = "R²", range = c(0, 1)),
        legend  = list(title = list(text = "SNP Type")),
        barmode = "relative"
      )
  })

  ## -- Download handler: exports all rows matching current filters --
  output$download_csv <- downloadHandler(
    filename = function() paste0("ld_partners_", query_snp(), ".csv"),
    content  = function(file) write.csv(partners_df(), file, row.names = FALSE)
  )

  ## -- Table: LD partners --
  output$ld_table <- renderDT({
    df    <- partners_df()
    req(nrow(df) > 0)

    acols <- ANCESTRY_COLS[[input$ancestry]]   # fix: reactive, not isolate()

    datatable(
      df,
      rownames   = FALSE,
      filter     = "top",
      selection  = "single",
      options    = list(
        pageLength = 20,
        scrollX    = TRUE,
        order      = list(list(1, "desc"))
      )
    ) %>%
      formatRound(c("r2", "A1FREQ", "BETA", "SE"), digits = 4) %>%
      formatRound("LOG10P", digits = 2) %>%
      formatStyle("Type",
        color      = styleEqual(
          c("VNTR region", "Non-repetitive region"),
          c(acols$vntr, acols$nonrep)       # fix: now reactive with ancestry
        ),
        fontWeight = styleEqual("VNTR region", "bold")
      ) %>%
      formatStyle("r2",
        background = styleColorBar(c(0, 1), "#dce8f5"),
        backgroundSize   = "98% 80%",
        backgroundRepeat = "no-repeat",
        backgroundPosition = "center"
      )
  })
}

shinyApp(ui, server)
