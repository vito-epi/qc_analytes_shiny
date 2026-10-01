required_packages <- c(
  "shiny",
  "dplyr",
  "purrr",
  "DT",
  "ggplot2"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  
  message(
    "Installing missing packages: ",
    paste(missing_packages, collapse = ", ")
  )
  
  install.packages(missing_packages)
}

invisible(
  lapply(required_packages, library, character.only = TRUE)
)

options(shiny.maxRequestSize = 500 * 1024^2)

required_objects <- c(
  "biomarker_lookup", "biomarker_overview", "variable_overview",
  "summary_statistics", "matrix_summary", "lod_loq_combinations",
  "imp_integrity", "meb_integrity", "plot_data", "distribution_data",
  "qq_data", "qq_reference", "flags"
)

imputation_mapping <- tibble(
  vartype = c("imp", "meb", "bin"),
  status_column = c("imp_status", "meb_status", "bin_status"),
  pct_column = c("imp_pct_complete", "meb_pct_complete", "bin_pct_complete")
)

correction_mapping <- tibble(
  vartype = c("imp_sg", "meb_sg", "imp_crt", "meb_crt", "imp_lip", "meb_lip"),
  status_column = c("imp_sg_status", "meb_sg_status", "imp_crt_status", "meb_crt_status", "imp_lip_status", "meb_lip_status"),
  pct_column = c("imp_sg_pct_complete", "meb_sg_pct_complete", "imp_crt_pct_complete", "meb_crt_pct_complete", "imp_lip_pct_complete", "meb_lip_pct_complete")
)

status_colour <- function(status, vartype = NA_character_) {
  if (identical(status, "PARTIAL") && vartype %in% c("meb", "bin")) return("red")
  case_when(
    status == "COMPLETE" ~ "darkgreen",
    status == "EXPECTED_EMPTY" ~ "black",
    status == "PARTIAL" ~ "#E69F00",
    status %in% c("UNEXPECTED_EMPTY", "VARIABLE_MISSING") ~ "red",
    TRUE ~ "#666666"
  )
}

status_label <- function(status, vartype = NA_character_) {
  if (identical(status, "PARTIAL") && vartype %in% c("meb", "bin")) return("partially complete; QC failure")
  case_when(
    status == "COMPLETE" ~ "complete",
    status == "EXPECTED_EMPTY" ~ "expected empty",
    status == "PARTIAL" ~ "partially complete",
    status == "UNEXPECTED_EMPTY" ~ "unexpectedly empty",
    status == "VARIABLE_MISSING" ~ "variable missing",
    TRUE ~ "not assessed"
  )
}

overall_status_colour <- function(status) {
  case_when(status == "PASS" ~ "darkgreen", status == "WARNING" ~ "#E69F00", status == "FAIL" ~ "red", TRUE ~ "#666666")
}

format_percentage <- function(x) {
  if (length(x) == 0 || is.na(x)) return("unavailable")
  paste0(format(x, trim = TRUE, nsmall = 1), "%")
}

make_qc_line <- function(varname, status, pct_complete, vartype) {
  percentage_text <- if (length(pct_complete) == 0 || is.na(pct_complete)) "completion unavailable" else paste0(format(pct_complete, trim = TRUE, nsmall = 1), "% complete")
  tags$li(style = paste0("color:", status_colour(status, vartype), ";font-weight:bold;"), paste0(varname, " (", percentage_text, "; ", status_label(status, vartype), ")"))
}

make_qc_lines <- function(expected_variables, biomarker_qc, mapping) {
  map(seq_len(nrow(mapping)), function(i) {
    current_type <- mapping$vartype[[i]]
    expected_row <- expected_variables %>% filter(vartype == current_type) %>% slice(1)
    if (nrow(expected_row) == 0) return(NULL)
    status_column <- mapping$status_column[[i]]
    pct_column <- mapping$pct_column[[i]]
    if (!status_column %in% names(biomarker_qc) || !pct_column %in% names(biomarker_qc)) return(NULL)
    make_qc_line(expected_row$Varname[[1]], biomarker_qc[[status_column]][[1]], biomarker_qc[[pct_column]][[1]], current_type)
  }) %>% compact()
}

make_progress_bar <- function(segments) {
  tags$div(
    style = "display:flex;height:18px;width:100%;border:1px solid #ccc;border-radius:3px;overflow:hidden;background:#f2f2f2;",
    lapply(segments, function(segment) tags$div(
      title = segment$title,
      style = paste0("width:", segment$width, "%;background:", segment$colour, ";min-width:", if (segment$width > 0) "2px" else "0", ";")
    ))
  )
}

ui <- fluidPage(
  tags$head(tags$style(HTML("
    .upload-card {max-width:760px;margin:40px auto;padding:24px;border:1px solid #ccc;border-radius:8px;background:#fafafa;}
    .sidebar-panel {max-height:calc(100vh - 80px);overflow-y:auto;}
    .compact-card {border:1px solid #ddd;border-radius:5px;padding:12px;margin-bottom:12px;background:#fafafa;}
    .compact-card h4,.compact-card h5 {margin-top:0;}
    .summary-stats-table table {font-size:12px;margin-bottom:0;}
    .summary-stats-table th,.summary-stats-table td {padding:5px !important;}
  "))),
  titlePanel("QC Dashboard of analyte data"),
  uiOutput("app_ui")
)

server <- function(input, output, session) {
  qc_objects <- reactiveVal(NULL)
  upload_error <- reactiveVal(NULL)

  observeEvent(input$qc_file, {
    req(input$qc_file)
    upload_error(NULL)
    tryCatch({
      x <- readRDS(input$qc_file$datapath)
      if (!is.list(x)) stop("The uploaded RDS must contain a named list.")
      missing_objects <- setdiff(required_objects, names(x))
      if (length(missing_objects) > 0) stop(paste("Missing QC objects:", paste(missing_objects, collapse = ", ")))
      qc_objects(x)
    }, error = function(e) {
      qc_objects(NULL)
      upload_error(conditionMessage(e))
    })
  })

  output$app_ui <- renderUI({
    if (is.null(qc_objects())) {
      return(div(
        class = "upload-card",
        h3("Open a QC object"),
        p("Select the qc_objects.rds file. The file is read only by this locally running R session."),
        fileInput("qc_file", "QC object (.rds)", accept = c(".rds", "application/octet-stream"), buttonLabel = "Browse...", placeholder = "No file selected"),
        if (!is.null(upload_error())) div(class = "alert alert-danger", tags$b("File could not be loaded: "), upload_error()),
        tags$hr(),
        p(style = "color:#666;font-size:12px;", "Close the Shiny app or browser tab when the review is finished.")
      ))
    }

    fluidRow(
      column(3,
        div(class = "sidebar-panel",
          actionButton("change_file", "Load another QC file", width = "100%"), br(), br(),
          selectizeInput("matrix", "Matrix", choices = NULL),
          uiOutput("matrix_summary_compact"),
          checkboxGroupInput("qc_status_filter", "QC status", choices = c("PASS", "WARNING", "FAIL"), selected = c("PASS", "WARNING", "FAIL"), inline = TRUE),
          uiOutput("filter_summary"),
          actionLink("select_problems", "Problems only"), tags$span(" | "), actionLink("select_all_statuses", "Show all"), br(), br(),
          selectizeInput("biomarker", "Biomarker", choices = NULL, options = list(placeholder = "Search by code or biomarker name...", maxOptions = 1000)),
          fluidRow(
            column(4, actionButton("previous_biomarker", "Previous", width = "100%")),
            column(4, div(style = "text-align:center;padding-top:7px;font-weight:bold;", textOutput("current_position", inline = TRUE))),
            column(4, actionButton("next_biomarker", "Next", width = "100%"))
          ), br(),
          div(class = "compact-card", h4("Biomarker"), uiOutput("biomarker_card"))
        )
      ),
      column(9,
             tabsetPanel(
             tabPanel("QC Summary", br(), fluidRow(
               
               column(6,
                      div(class = "compact-card",
                          uiOutput("summary_table")
                      ),
                      
                      div(class = "compact-card summary-stats-table",
                          h4("Summary statistics"),
                          h5("Continuous variables"),
                          tableOutput("continuous_summary_table"),
                          tags$hr(style = "margin:12px 0;"),
                          h5("Detection variable"),
                          tableOutput("bin_summary_table")
                      ),
                      
                      div(class = "compact-card",
                          uiOutput("qc_info")
                      )
               ),
               
               column(6,
                      div(class = "compact-card",
                          h4("LOD / LOQ"),
                          uiOutput("lod_loq_message"),
                          tableOutput("lod_loq_table")
                      ),
                      
                      div(class = "compact-card",
                          uiOutput("integrity_overview")
                      ),
                      
                      # div(class = "compact-card",
                      #     h4("Distribution QC"),
                      #     plotOutput("qq_plot_summary",
                      #                height = "300px")
                      # )
                      div(
                        class = "compact-card",
                        
                        h4("Distribution QC"),
                        
                        fluidRow(
                          
                          column(
                            6,
                            plotOutput(
                              "distribution_plot_summary",
                              height = "280px"
                            )
                          ),
                          
                          column(
                            6,
                            plotOutput(
                              "qq_plot_summary",
                              height = "280px"
                            )
                          )
                          
                        )
                      )
                      
               )
               
             )),
             # tabPanel(
             #   "Plots",
             #   br(),
             #   
             #   fluidRow(
             #     column(
             #       12,
             #       plotOutput(
             #         "distribution_plot",
             #         height = "420px"
             #       )
             #     )
             #   ),
             #   
             #   fluidRow(
             #     column(
             #       12,
             #       plotOutput(
             #         "qq_plot",
             #         height = "420px"
             #       )
             #     )
             #   )
             # ),
             tabPanel("All biomarkers", br(), DTOutput("overview_table"))
        )
      )
    )
  })

  observeEvent(input$change_file, { qc_objects(NULL); upload_error(NULL) })

  obj <- function(name) {
    req(qc_objects())
    qc_objects()[[name]]
  }

  observeEvent(qc_objects(), {
    req(qc_objects())
    matrices <- sort(unique(obj("biomarker_overview")$matrix))
    updateSelectizeInput(session, "matrix", choices = matrices, selected = matrices[[1]], server = FALSE)
  })

  observeEvent(input$select_problems, updateCheckboxGroupInput(session, "qc_status_filter", selected = c("WARNING", "FAIL")))
  observeEvent(input$select_all_statuses, updateCheckboxGroupInput(session, "qc_status_filter", selected = c("PASS", "WARNING", "FAIL")))

  matrix_biomarkers <- reactive({ req(input$matrix); obj("biomarker_overview") %>% filter(matrix == input$matrix) })
  filtered_biomarkers <- reactive({
    statuses <- input$qc_status_filter
    validate(need(length(statuses) > 0, "Select at least one QC status."))
    matrix_biomarkers() %>% filter(Overall_QC_Status %in% statuses)
  })
  biomarker_choices <- reactive({
    req(input$matrix)
    statuses <- input$qc_status_filter
    validate(need(length(statuses) > 0, "Select at least one QC status."))
    allowed <- obj("biomarker_overview") %>% filter(matrix == input$matrix, Overall_QC_Status %in% statuses) %>% pull(BaseVarname)
    obj("biomarker_lookup") %>%
      filter(matrix == input$matrix, BaseVarname %in% allowed) %>%
      mutate(
        Description = if_else(is.na(Description) | Description == "", BaseVarname, Description),
        Unit = if_else(is.na(Unit) | Unit == "", "unit unavailable", Unit),
        display_label = paste0(Description, " [", Unit, "] (", BaseVarname, ")")
      ) %>% arrange(Tabname, Order, Suborder, Description, BaseVarname) %>% distinct(BaseVarname, .keep_all = TRUE)
  })
  available_biomarkers <- reactive(biomarker_choices()$BaseVarname)

  observeEvent(list(input$matrix, input$qc_status_filter), {
    req(qc_objects(), input$matrix)
    choices_df <- biomarker_choices()
    if (nrow(choices_df) == 0) {
      updateSelectizeInput(session, "biomarker", choices = character(), selected = character(), server = TRUE)
      return()
    }
    choices <- setNames(choices_df$BaseVarname, choices_df$display_label)
    current <- isolate(input$biomarker)
    selected <- if (!is.null(current) && current %in% choices_df$BaseVarname) current else choices_df$BaseVarname[[1]]
    updateSelectizeInput(session, "biomarker", choices = choices, selected = selected, server = FALSE)
  }, ignoreInit = TRUE)

  observeEvent(input$previous_biomarker, {
    biomarkers <- available_biomarkers(); req(length(biomarkers) > 0)
    i <- match(input$biomarker, biomarkers); if (is.na(i)) i <- 1
    i <- if (i == 1) length(biomarkers) else i - 1
    updateSelectizeInput(session, "biomarker", selected = biomarkers[[i]], server = FALSE)
  })
  observeEvent(input$next_biomarker, {
    biomarkers <- available_biomarkers(); req(length(biomarkers) > 0)
    i <- match(input$biomarker, biomarkers); if (is.na(i)) i <- 0
    i <- if (i >= length(biomarkers)) 1 else i + 1
    updateSelectizeInput(session, "biomarker", selected = biomarkers[[i]], server = FALSE)
  })

  output$current_position <- renderText({
    biomarkers <- available_biomarkers(); if (length(biomarkers) == 0) return("0 / 0")
    i <- match(input$biomarker, biomarkers); if (is.na(i)) i <- 0
    paste0(i, " / ", length(biomarkers))
  })

  current_qc <- reactive({
    req(input$matrix, input$biomarker)
    out <- obj("biomarker_overview") %>% filter(matrix == input$matrix, BaseVarname == input$biomarker)
    validate(need(nrow(out) == 1, "No unique QC record found for the selected biomarker.")); out
  })
  current_imp_integrity <- reactive({
    out <- obj("imp_integrity") %>% filter(matrix == input$matrix, BaseVarname == input$biomarker)
    validate(need(nrow(out) == 1, "No IMP integrity information available.")); out
  })
  current_meb_integrity <- reactive({
    out <- obj("meb_integrity") %>% filter(matrix == input$matrix, BaseVarname == input$biomarker)
    validate(need(nrow(out) == 1, "No MEB integrity information available.")); out
  })
  current_summary_statistics <- reactive({
    obj("summary_statistics") %>% filter(matrix == input$matrix, BaseVarname == input$biomarker, vartype %in% c("raw", "imp", "meb", "bin")) %>% mutate(vartype = factor(vartype, levels = c("raw", "imp", "meb", "bin"))) %>% arrange(vartype)
  })
  selected_variables <- reactive({
    type_order <- c("raw", "lod", "loq", "imp", "meb", "bin", "imp_sg", "meb_sg", "imp_crt", "meb_crt", "imp_lip", "meb_lip", "other")
    obj("variable_overview") %>% filter(matrix == input$matrix, BaseVarname == input$biomarker) %>% mutate(vartype_order = match(vartype, type_order)) %>% arrange(vartype_order, Varname) %>% select(-vartype_order)
  })
  current_lod_loq <- reactive(obj("lod_loq_combinations") %>% filter(matrix == input$matrix, BaseVarname == input$biomarker) %>% arrange(desc(N), LOD, LOQ))
  current_plot_data <- reactive(obj("plot_data") %>% filter(matrix == input$matrix, BaseVarname == input$biomarker))
  current_distribution_data <- reactive(obj("distribution_data") %>% filter(matrix == input$matrix, BaseVarname == input$biomarker))
  current_qq_data <- reactive(obj("qq_data") %>% filter(matrix == input$matrix, BaseVarname == input$biomarker))
  current_qq_reference <- reactive(obj("qq_reference") %>% filter(matrix == input$matrix, BaseVarname == input$biomarker))

  output$matrix_summary_compact <- renderUI({
    s <- obj("matrix_summary") %>% filter(matrix == input$matrix)
    validate(need(nrow(s) == 1, "No matrix summary available."))
    tags$p(style = "font-size:12px;line-height:1.35;", paste0(s$N_records[[1]], " records | ", s$N_biomarkers[[1]], " biomarkers"), tags$br(), tags$span(style = "color:red;font-weight:bold;", paste0(s$N_QC_fail[[1]], " FAIL")), " | ", tags$span(style = "color:#E69F00;font-weight:bold;", paste0(s$N_QC_warning[[1]], " WARNING")), " | ", tags$span(style = "color:darkgreen;font-weight:bold;", paste0(s$N_QC_pass[[1]], " PASS")))
  })
  output$filter_summary <- renderUI({
    statuses <- input$qc_status_filter
    if (length(statuses) == 0) return(tags$p(style = "color:red;font-weight:bold;", "No QC status selected."))
    tags$p(style = "color:#666;font-size:12px;", paste0(nrow(filtered_biomarkers()), " of ", nrow(matrix_biomarkers()), " biomarkers shown"))
  })
  output$biomarker_card <- renderUI({
    s <- current_qc()
    tags$div(tags$b(s$Description[[1]]), tags$br(), tags$span(style = "color:#666;", paste0(s$BaseVarname[[1]], " | ", s$Unit[[1]])), tags$br(), paste0(s$N_non_missing[[1]], " / ", s$N_total[[1]], " non-missing"), tags$br(), tags$span(style = paste0("color:", overall_status_colour(s$Overall_QC_Status[[1]]), ";font-weight:bold;"), s$Overall_QC_Status[[1]]))
  })
  output$summary_table <- renderUI({
    s <- current_qc(); pct_ok <- isTRUE(s$Pct_detected[[1]] >= 30); unique_ok <- isTRUE(s$N_unique_detected[[1]] >= 10); imp_expected <- isTRUE(s$IMP_expected[[1]])
    tagList(h4("Imputation criteria"),tags$p(
      tags$b("N values above LOD/LOQ: "),
      paste0(
        s$N_detected[[1]],
        " of ",
        s$N_non_missing[[1]],
        " records"
      )
    ), tags$p(tags$b("% above LOD/LOQ (>=30%): "), tags$span(style = if (pct_ok) "color:darkgreen;font-weight:bold;" else "color:red;font-weight:bold;", format_percentage(s$Pct_detected[[1]]))), tags$p(tags$b("Unique values above LOD/LOQ (>=10): "), tags$span(style = if (unique_ok) "color:darkgreen;font-weight:bold;" else "color:red;font-weight:bold;", s$N_unique_detected[[1]])), tags$p(tags$b("IMP expected: "), ifelse(imp_expected, "YES", "NO")), tags$p(tags$b("Overall QC status: "), tags$span(style = paste0("color:", overall_status_colour(s$Overall_QC_Status[[1]]), ";font-weight:bold;"), s$Overall_QC_Status[[1]])))
  })
  output$continuous_summary_table <- renderTable({
    out <- current_summary_statistics() %>% filter(vartype %in% c("raw", "imp", "meb")); validate(need(nrow(out) > 0, "No continuous summary statistics available."))
    out %>% mutate(Variable = as.character(vartype), across(any_of(c("Min", "P25", "Median", "Mean", "P75", "Max")), ~ round(.x, 3))) %>% select(Variable, N, Missing = N_missing, Min, P25, Median, Mean, P75, Max)
  }, striped = TRUE, bordered = TRUE, spacing = "xs", width = "100%", na = "")
  output$bin_summary_table <- renderTable({
    out <- current_summary_statistics() %>% filter(vartype == "bin"); validate(need(nrow(out) > 0, "No detection variable available."))
    out %>% mutate(Variable = as.character(vartype), Pct_0 = if_else(N_0 + N_1 > 0, round(100 * N_0 / (N_0 + N_1), 1), NA_real_), Pct_1 = round(Pct_1, 1)) %>% select(Variable, N, Missing = N_missing, `N = 0` = N_0, `% = 0` = Pct_0, `N = 1` = N_1, `% = 1` = Pct_1)
  }, striped = TRUE, bordered = TRUE, spacing = "xs", width = "100%", na = "")
  output$qc_info <- renderUI({
    s <- current_qc(); ev <- selected_variables(); imp_lines <- make_qc_lines(ev, s, imputation_mapping); corr_lines <- make_qc_lines(ev, s, correction_mapping)
    tagList(h4("QC summary"), 
            # tags$p(tags$b("Imputed observations: "), paste0(s$N_imputed[[1]], " (", format_percentage(s$Pct_imputed[[1]]), ")")), 
            tags$p(tags$b("Censored observations: "), paste0(s$N_censored[[1]], " (", format_percentage(s$Pct_non_detected[[1]]), ")")), h5("Imputation variables"), if (length(imp_lines) > 0) tags$ul(imp_lines) else tags$p("No imputation variables defined."), h5("Correction variables"), if (length(corr_lines) > 0) tags$ul(corr_lines) else tags$p("No correction variables defined."))
  })
  output$lod_loq_message <- renderUI({
    s <- current_qc(); out <- current_lod_loq()
    if (nrow(out) == 0) return(tags$p(style = "color:#666;", "No LOD/LOQ information available."))
    tagList(tags$p(paste0(s$N_lod_loq_combinations[[1]], " unique combination", ifelse(s$N_lod_loq_combinations[[1]] == 1, "", "s"))), if (isTRUE(s$Warning_missing_LOD_LOQ[[1]])) tags$p(style = "color:red;font-weight:bold;", paste0("Warning: ", s$N_missing_LOD_LOQ[[1]], " observation(s) have a raw value but both LOD and LOQ are missing.")))
  })
  output$lod_loq_table <- renderTable({
    out <- current_lod_loq(); validate(need(nrow(out) > 0, "")); out %>% select(LOD, LOQ, N, Pct) %>% rename(`Percentage (%)` = Pct)
  }, striped = TRUE, bordered = TRUE, spacing = "xs", width = "100%")
  output$integrity_overview <- renderUI({
    imp <- current_imp_integrity(); meb <- current_meb_integrity()
    measured_total <- imp$IMP_measured_total[[1]]; censored_total <- imp$IMP_censored_total[[1]]; raw_total <- measured_total + censored_total
    safe_count <- function(x) if (length(x) == 0 || is.na(x)) 0 else x
    pct_width <- function(x) if (is.na(raw_total) || raw_total <= 0 || is.na(x)) 0 else 100 * x / raw_total
    imp_expected <- isTRUE(imp$IMP_expected[[1]]); imp_status <- imp$IMP_integrity_status[[1]]
    ir <- safe_count(imp$IMP_measured_retained[[1]]); ii <- safe_count(imp$IMP_censored_imputed[[1]]); im <- safe_count(imp$IMP_measured_missing[[1]]) + safe_count(imp$IMP_censored_missing[[1]])
    mr <- safe_count(meb$MEB_measured_retained[[1]]); ms <- safe_count(meb$MEB_censored_substituted[[1]]); mm <- safe_count(meb$MEB_measured_missing[[1]]) + safe_count(meb$MEB_censored_missing[[1]])
    tagList(
      h4("IMP and MEB integrity"),
      tags$div(style = "display:flex;justify-content:space-between;", tags$b("Raw"), tags$span(paste0(measured_total, " measured | ", censored_total, " censored"))),
      make_progress_bar(list(list(width = pct_width(measured_total), colour = "#2C7FB8", title = paste0(measured_total, " measured")), list(width = pct_width(censored_total), colour = "#D95F0E", title = paste0(censored_total, " censored")))),
      tags$hr(style = "margin:12px 0 10px 0;"),
      tags$div(style = "display:flex;justify-content:space-between;", tags$b("IMP"), tags$span(if (!imp_expected || identical(imp_status, "NOT_EXPECTED")) "not expected" else paste0(ir, " of measured | ", ii, " of censored"))),
      if (!imp_expected || identical(imp_status, "NOT_EXPECTED")) tags$p(style = "color:#666;font-size:12px;", "The _imp variable is not assessed because imputation is not expected.") else tagList(make_progress_bar(list(list(width = pct_width(ir), colour = "#31A354", title = paste0(ir, " retained")), list(width = pct_width(ii), colour = "#756BB1", title = paste0(ii, " imputed")), list(width = pct_width(im), colour = "#CB181D", title = paste0(im, " missing from _imp")))), tags$p(style = paste0("font-size:12px;font-weight:bold;color:", if (im > 0) "red" else "darkgreen", ";"), if (im > 0) paste0(im, " observation(s) missing from _imp.") else "all measured values have a non-NA value; all censored values have a non-NA value")),
      tags$div(style = "height:12px;"),
      tags$div(style = "display:flex;justify-content:space-between;", tags$b("MEB"), tags$span(paste0(mr, " of measured | ", ms, " of censored"))),
      make_progress_bar(list(list(width = pct_width(mr), colour = "#31A354", title = paste0(mr, " retained")), list(width = pct_width(ms), colour = "#756BB1", title = paste0(ms, " substituted")), list(width = pct_width(mm), colour = "#CB181D", title = paste0(mm, " missing from _meb")))),
      tags$p(style = paste0("font-size:12px;font-weight:bold;color:", if (mm > 0) "red" else "darkgreen", ";"), if (mm > 0) paste0(mm, " observation(s) missing from _meb.") else "all measured values have a non-NA value; all censored values have a non-NA value")
    )
  })
  output$overview_table <- renderDT({
    out <- filtered_biomarkers() %>% arrange(factor(Overall_QC_Status, levels = c("FAIL", "WARNING", "PASS")), BaseVarname) %>% select(BaseVarname, Description, Unit, N_total, N_non_missing, N_detected, Pct_detected, N_unique_detected, IMP_expected, N_censored, N_imputed, Pct_imputed, N_lod_loq_combinations, Warning_missing_LOD_LOQ, imp_status, IMP_integrity_status, meb_status, MEB_integrity_status, bin_status, Overall_QC_Status)
    datatable(out, rownames = FALSE, filter = "top", options = list(pageLength = 25, scrollX = TRUE)) %>% formatStyle("Overall_QC_Status", color = styleEqual(c("PASS", "WARNING", "FAIL"), c("darkgreen", "#E69F00", "red")), fontWeight = "bold")
  })
  # output$hist_plot <- renderPlot({
  #   d <- current_plot_data(); validate(need(nrow(d) > 0, "No positive measured or imputed values available."))
  #   ggplot(d, aes(x = value, fill = source)) + geom_histogram(alpha = 0.5, bins = 30, position = "identity") + scale_x_log10() + scale_fill_manual(values = c(Measured = "#2C7FB8", Imputed = "#D95F0E"), drop = FALSE) + theme_minimal() + labs(title = paste0(input$biomarker, ": measured versus imputed values"), x = "Value (log scale)", y = "Count", fill = NULL)
  # })
  # output$density_plot <- renderPlot({
  #   d <- current_distribution_data(); validate(need(nrow(d) >= 10, paste0("Too few positive non-missing observations in ", input$biomarker, "_imp.")))
  #   ggplot(d, aes(x = value)) + geom_density(fill = "#8E63CE", colour = "#8E63CE", alpha = 0.3, linewidth = 1) + scale_x_log10() + theme_minimal() + labs(title = paste0("Full distribution of ", input$biomarker, "_imp"), x = paste0(input$biomarker, "_imp (log scale)"), y = "Density")
  # })
  output$distribution_plot <- renderPlot({
    
    grouped_data <- current_plot_data()
    full_data <- current_distribution_data()
    
    validate(
      need(
        nrow(grouped_data) > 0,
        "No positive measured or imputed values available."
      ),
      need(
        nrow(full_data) >= 10,
        paste0(
          "Too few positive non-missing observations in ",
          input$biomarker,
          "_imp."
        )
      )
    )
    
    grouped_data <- grouped_data %>%
      filter(
        is.finite(value),
        value > 0,
        source %in% c("Measured", "Imputed")
      ) %>%
      mutate(
        log_value = log10(value),
        source = factor(
          source,
          levels = c("Measured", "Imputed")
        )
      )
    
    full_data <- full_data %>%
      filter(
        is.finite(value),
        value > 0
      ) %>%
      mutate(
        log_value = log10(value)
      )
    
    validate(
      need(
        nrow(grouped_data) > 0,
        "No positive measured or imputed values available."
      ),
      need(
        nrow(full_data) >= 10,
        "Too few positive non-missing values for the full distribution."
      ),
      need(
        dplyr::n_distinct(full_data$log_value) >= 2,
        "At least two distinct positive values are needed for a density curve."
      )
    )
    
    number_of_bins <- 30
    
    log_range <- range(
      c(
        grouped_data$log_value,
        full_data$log_value
      ),
      finite = TRUE
    )
    
    bin_width <- diff(log_range) / number_of_bins
    
    validate(
      need(
        is.finite(bin_width) && bin_width > 0,
        "Insufficient variation to construct a distribution plot."
      )
    )
    
    density_estimate <- density(
      full_data$log_value,
      from = log_range[[1]],
      to = log_range[[2]],
      n = 512,
      na.rm = TRUE
    )
    
    density_data <- tibble(
      log_value = density_estimate$x,
      count = density_estimate$y *
        nrow(full_data) *
        bin_width
    )
    
    exponent_min <- floor(log_range[[1]])
    exponent_max <- ceiling(log_range[[2]])
    
    log_breaks <- seq(
      exponent_min,
      exponent_max,
      by = 1
    )
    
    log_labels <- format(
      10^log_breaks,
      scientific = FALSE,
      trim = TRUE,
      big.mark = ","
    )
    
    ggplot() +
      
      geom_histogram(
        data = grouped_data,
        aes(
          x = log_value,
          fill = source
        ),
        bins = number_of_bins,
        position = "identity",
        alpha = 0.60,
        colour = NA
      ) +
      
      geom_line(
        data = density_data,
        aes(
          x = log_value,
          y = count,
          colour = "Full distribution"
        ),
        linewidth = 1.2
      ) +
      
      scale_fill_manual(
        values = c(
          Measured = "#2C7FB8",
          Imputed = "#D95F0E"
        ),
        drop = FALSE
      ) +
      
      scale_colour_manual(
        values = c(
          `Full distribution` = "#8E63CE"
        )
      ) +
      
      scale_x_continuous(
        breaks = log_breaks,
        labels = log_labels,
        expand = expansion(mult = c(0.02, 0.03))
      ) +
      
      scale_y_continuous(
        expand = expansion(mult = c(0, 0.05))
      ) +
      
      guides(
        fill = guide_legend(
          title = NULL,
          order = 1,
          override.aes = list(alpha = 0.60)
        ),
        colour = guide_legend(
          title = NULL,
          order = 2,
          override.aes = list(linewidth = 1.2)
        )
      ) +
      
      theme_minimal() +
      
      theme(
        legend.position = "right",
        panel.grid.minor = element_blank()
      ) +
      
      labs(
        title = paste0(
          input$biomarker,
          ": measured, imputed and full distribution"
        ),
        subtitle = paste0(
          "Bars show measured and imputed observations; ",
          "the purple line shows the complete ",
          input$biomarker,
          "_imp distribution."
        ),
        x = "Value (log10 scale)",
        y = "Count",
        fill = NULL,
        colour = NULL
      )
  })
  
  
  output$distribution_plot_summary <- renderPlot({
    
    grouped_data <- current_plot_data()
    full_data <- current_distribution_data()
    
    validate(
      need(
        nrow(grouped_data) > 0,
        "No positive measured or imputed values available."
      ),
      need(
        nrow(full_data) >= 10,
        "Too few positive non-missing observations."
      )
    )
    
    grouped_data <- grouped_data %>%
      filter(
        is.finite(value),
        value > 0,
        source %in% c("Measured", "Imputed")
      ) %>%
      mutate(
        log_value = log10(value),
        source = factor(
          source,
          levels = c("Measured", "Imputed")
        )
      )
    
    full_data <- full_data %>%
      filter(
        is.finite(value),
        value > 0
      ) %>%
      mutate(
        log_value = log10(value)
      )
    
    validate(
      need(
        nrow(grouped_data) > 0,
        "No positive measured or imputed values available."
      ),
      need(
        nrow(full_data) >= 10,
        "Too few positive non-missing values for the full distribution."
      ),
      need(
        dplyr::n_distinct(full_data$log_value) >= 2,
        "At least two distinct positive values are needed."
      )
    )
    
    number_of_bins <- 30
    
    log_range <- range(
      c(
        grouped_data$log_value,
        full_data$log_value
      ),
      finite = TRUE
    )
    
    bin_width <- diff(log_range) / number_of_bins
    
    validate(
      need(
        is.finite(bin_width) && bin_width > 0,
        "Insufficient variation to construct a distribution plot."
      )
    )
    
    density_estimate <- density(
      full_data$log_value,
      from = log_range[[1]],
      to = log_range[[2]],
      n = 512,
      na.rm = TRUE
    )
    
    density_data <- tibble(
      log_value = density_estimate$x,
      count = density_estimate$y *
        nrow(full_data) *
        bin_width
    )
    
    exponent_min <- floor(log_range[[1]])
    exponent_max <- ceiling(log_range[[2]])
    
    log_breaks <- seq(
      exponent_min,
      exponent_max,
      by = 1
    )
    
    log_labels <- format(
      10^log_breaks,
      scientific = FALSE,
      trim = TRUE,
      big.mark = ","
    )
    
    ggplot() +
      
      geom_histogram(
        data = grouped_data,
        aes(
          x = log_value,
          fill = source
        ),
        bins = number_of_bins,
        position = "identity",
        alpha = 0.60,
        colour = NA
      ) +
      
      geom_line(
        data = density_data,
        aes(
          x = log_value,
          y = count,
          colour = "Full distribution"
        ),
        linewidth = 1
      ) +
      
      scale_fill_manual(
        values = c(
          Measured = "#2C7FB8",
          Imputed = "#D95F0E"
        ),
        drop = FALSE
      ) +
      
      scale_colour_manual(
        values = c(
          `Full distribution` = "#8E63CE"
        )
      ) +
      
      scale_x_continuous(
        breaks = log_breaks,
        labels = log_labels,
        expand = expansion(mult = c(0.02, 0.03))
      ) +
      
      scale_y_continuous(
        expand = expansion(mult = c(0, 0.05))
      ) +
      
      guides(
        fill = guide_legend(
          title = NULL,
          order = 1,
          override.aes = list(alpha = 0.60)
        ),
        colour = guide_legend(
          title = NULL,
          order = 2,
          override.aes = list(linewidth = 1)
        )
      ) +
      
      theme_minimal(base_size = 10) +
      
      theme(
        legend.position = "bottom",
        legend.box = "vertical",
        legend.spacing.y = unit(0, "pt"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
        panel.grid.minor = element_blank(),
        plot.margin = margin(5, 5, 5, 5)
      ) +
      
      labs(
        title = "Measured, imputed and full distribution",
        x = "Value (log10 scale)",
        y = "Count",
        fill = NULL,
        colour = NULL
      )
  })
  
  output$qq_plot <- renderPlot({
    d <- current_qq_data(); ref <- current_qq_reference()
    validate(need(nrow(d) >= 10, "Too few positive non-missing observations."), need(nrow(ref) == 1, "No Q-Q reference line available."))
    ggplot(d, aes(x = theoretical, y = log_value, colour = source)) + geom_point(size = 2, alpha = 0.8) + geom_abline(slope = ref$qq_slope[[1]], intercept = ref$qq_intercept[[1]], colour = "black", linewidth = 1) + scale_colour_manual(values = c(Measured = "#2C7FB8", Imputed = "#D95F0E"), drop = FALSE) + theme_minimal() + labs(title = paste0("Normal Q-Q assessment of log(", input$biomarker, "_imp)"), x = "Theoretical normal quantiles", y = "Sample quantiles", colour = NULL)
  })
  output$qq_plot_summary <- renderPlot({
    
    d <- current_qq_data()
    ref <- current_qq_reference()
    
    validate(
      need(
        nrow(d) >= 10,
        "Too few positive non-missing observations."
      ),
      need(
        nrow(ref) == 1,
        "No Q-Q reference line available."
      )
    )
    
    ggplot(
      d,
      aes(
        x = theoretical,
        y = log_value,
        colour = source
      )
    ) +
      geom_point(
        size = 2,
        alpha = 0.8
      ) +
      geom_abline(
        slope = ref$qq_slope[[1]],
        intercept = ref$qq_intercept[[1]],
        colour = "black",
        linewidth = 1
      ) +
      scale_colour_manual(
        values = c(
          Measured = "#2C7FB8",
          Imputed = "#D95F0E"
        ),
        drop = FALSE
      ) +
      theme_minimal() +
      labs(
        title = paste0(
          "Normal Q-Q assessment of log(",
          input$biomarker,
          "_imp)"
        ),
        x = "Theoretical normal quantiles",
        y = "Sample quantiles",
        colour = NULL
      )
    
  })
  
}

shinyApp(ui = ui, server = server)
