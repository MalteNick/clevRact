.build_ui <- function() {
  fluidPage(
    useShinyjs(),
    tags$head(
      tags$style(HTML("
        body {
          overflow-y: scroll;
        }
        ::-webkit-scrollbar {
          width: 0px;
          background: transparent;
        }
        html {
          scrollbar-width: none;
        }
    ")),
      tags$link(
        rel = "icon",
        type = "image/png",
        href = "clevRact-assets/clevRact.png"
      )
    ),
    title = "clevRact",
    theme = shinytheme("flatly"),
    titlePanel(div(
      img(
        height = 100,
        width = 90,
        src = "clevRact-assets/clevRact.png",
        class = "pull-left",
        style = "margin-right: 10px;"
      ),
      "clevRact"
    )),
    tabsetPanel(
      tabPanel(
        "Input Data",
        sidebarLayout(
          sidebarPanel(
            width = 3,
            div(
              style = "text-align: justify;",
              div(style = "font-size: 25px;", strong("Input Data Editor")),
              paste(
                "This tab is for showcasing and editing the raw data input",
                "used in clevRact. You can either create an all-new dataset",
                "using the"
              ),
              strong("New Dataset"),
              paste(
                "button, or load an existing dataset by providing a .csv",
                "file via the"
              ),
              strong("Load Data"),
              paste(
                "button. The .csv files require at minimum the following",
                "columns, using a semicolon as separator:"
              ),
              tags$ul(
                tags$li(
                  strong("Patient ID"),
                  " (Any nonempty unique char without semicolon)"
                ),
                tags$li(
                  strong("Variant"),
                  "  (Any nonempty char without semicolon)"
                ),
                tags$li(strong("Clone"), " (Integer ID other than 0)"),
                tags$li(strong("Parent"), " (Integer ID with 0 for germline)"),
                tags$li(strong("Chromosome"), " (1-22, X, Y, NA)"),
                tags$li(strong("Start"), " (Valid numeric bp position or NA)"),
                tags$li(strong("End"), " (Valid numeric bp position or NA)"),
                tags$li(
                  strong("Variant Type"),
                  " (CNV_del, CNV_dup, CNV_loh, SNV/INDEL or NA)"
                )
              ),
              paste(
                "If data on different time points is to be provided, it",
                "needs to be given as TRUE or FALSE in columns named"
              ),
              strong("Ti"),
              " where ",
              strong("i"),
              " the timepoint the data are taken from, counting from 1."
            ), br(),
            actionButton(
              "change_timepoints",
              "Change Number of Timepoints",
              class = "sidebar-btn"
            ),
            div(
              style = "display: flex; gap: 10px; margin-top: 10px;",
              actionButton(
                "add_patient",
                "Add Patient",
                class = "sidebar-btn"
              ),
              actionButton(
                "delete_patient",
                "Delete Patient",
                class = "sidebar-btn"
              )
            ),
            div(
              style = "display: flex; gap: 10px;",
              actionButton(
                "add_variant",
                "Add Variant",
                class = "sidebar-btn"
              ),
              actionButton(
                "delete_selection",
                "Delete Variants",
                class = "sidebar-btn"
              )
            ),
            actionButton(
              "enable_disable_selection",
              "Toggle Variant Selection",
              class = "sidebar-btn"
            ),
          ),
          mainPanel(
            height = 10, width = 9,
            fluidRow(
              style = paste(
                "margin-bottom: 5px; background-color: #ecf0f1;",
                "padding-top: 5px; padding-right: 15px"
              ),
              column(
                2,
                div(
                  style = "text-align: center; padding-top: 20px;",
                  actionButton(
                    "new_dataset",
                    "New Dataset",
                    class = "sidebar-btn",
                    width = "100%"
                  )
                )
              ),
              column(
                2,
                div(
                  style = "text-align: center; padding-top: 20px;",
                  actionButton(
                    "load_dataset",
                    "Load Dataset",
                    class = "sidebar-btn",
                    width = "100%"
                  )
                )
              ),
              column(
                2,
                div(
                  style = "text-align: center; padding-top: 20px;",
                  downloadButton("downloadData", "Save Dataset",
                    style = paste(
                      "display: block; width: 100%; margin-bottom: 10px;",
                      "height: 40px; text-align: center; font-size: 14px;",
                      "padding-top: 8px; box-sizing: border-box;"
                    )
                  )
                )
              ),
              column(
                2,
                div(
                  style = "padding-top: 20px;",
                  fluidRow(
                    column(
                      6,
                      actionButton("previous_patient", NULL,
                        icon = icon("chevron-up"),
                        class = "sidebar-btn",
                        style = "width: 100%; height: 40px;"
                      )
                    ),
                    column(
                      6,
                      actionButton("next_patient", NULL,
                        icon = icon("chevron-down"),
                        class = "sidebar-btn",
                        style = "width: 100%; height: 40px;"
                      )
                    )
                  )
                )
              ),
              column(
                2,
                uiOutput("id_selector")
              ),
              column(
                2,
                div(
                  style = paste(
                    "padding-top: 23px; width: 200px;",
                    "padding-right: 0px;"
                  ),
                  shinyWidgets::switchInput(
                    inputId = "toggle_view",
                    label = "Graph View",
                    onLabel = "On",
                    offLabel = "Off",
                    labelWidth = "100%",
                    width = "100%",
                    size = "normal",
                    inline = TRUE
                  )
                )
              )
            ),
            uiOutput("input_data_main_panel"),
            style = "padding: 0;"
          )
        )
      ),
      tabPanel(
        "Processed Data",
        sidebarLayout(
          sidebarPanel(
            width = 3,
            div(
              style = "text-align: justify;",
              div(
                style = "font-size: 25px;",
                strong("Processed Data Overview")
              ),
              paste(
                "This tab is for showcasing the preprocessed dataset for",
                "analysis before further Signature and Trajectory",
                "computation. It is generated from the input data which can",
                "be changed in the Input Data tab and the settings below.",
                "You can edit the data entries or the settings and press",
                "the Generate Processed Data button for (re)processing.",
                "Editing the input data without reprocessing will have no",
                "effect on the processed data or the computed statistics",
                "and signatures."
              )
            ), br(),
            div(style = "font-size: 25px;", strong("Settings")),
            selectInput(
              inputId = "dropdown_variant_generation",
              label = "Generate variant labels?",
              choices = c("No", "Yes, using GRCh37", "Yes, using GRCh38"),
              selected = "Yes, using GRCh37"
            ),
            selectInput(
              inputId = "dropdown_summary_snvs",
              label = "- Information detail for SNVs/INDEL annotations:",
              choices = c("Summarize on gene level", "Use exact SNVs"),
              selected = "Summarize on gene level"
            ),
            selectInput(
              inputId = "dropdown_summary_cnvs",
              label = "- Information detail for CNV annotations:",
              choices = c(
                "Summarize on chromosome level",
                "Chromosome, p/q arm, in p/q arm",
                "Use exact CNVs"
              ),
              selected = "Summarize on chromosome level"
            ),
            selectInput(
              inputId = "dropdown_handler_reoccurrence",
              label = paste(
                "Handler for reoccurrence of identical variant labels",
                "within sequence:"
              ),
              choices = c(
                "No differentiation",
                "Enumeration",
                "Enumeration for SNVs/INDELs only",
                "Enumeration for CNVs only"
              ),
              selected = "No differentiation"
            ),
            radioButtons(
              inputId = "dropdown_cnv_clone_summary",
              label = paste(
                "Summarize same-type CNVs within same clone and chromosome?"
              ),
              choices = c("Yes" = TRUE, "No" = FALSE),
              selected = TRUE,
              inline = TRUE
            )
          ),
          mainPanel(
            height = 10,
            width = 9,
            fluidRow(
              style = paste(
                "margin-bottom: 5px; background-color: #ecf0f1;",
                "padding-top: 5px; padding-right: 15px"
              ),
              column(
                3,
                div(
                  style = "text-align: center; padding-top: 15px;",
                  actionButton(
                    "preprocess_data",
                    "Generate Processed Data",
                    class = "btn btn-success",
                    width = "100%"
                  )
                )
              ),
              column(
                3,
                uiOutput("processed_id_selector")
              ),
              column(
                3,
                div(
                  style = "padding-top: 20px;",
                  actionButton(
                    "previous_processed_patient",
                    "Previous",
                    class = "sidebar-btn",
                    style = "width: 100%;"
                  )
                )
              ),
              column(
                3,
                div(
                  style = "padding-top: 20px;",
                  actionButton(
                    "next_processed_patient",
                    "Next",
                    class = "sidebar-btn",
                    style = "width: 100%;"
                  )
                )
              )
            ),
            uiOutput("processed_data_main_panel"),
            style = "padding: 0;"
          )
        )
      ),
      tabPanel(
        "Data Statistics",
        sidebarLayout(
          sidebarPanel(
            width = 3,
            div(style = "font-size: 25px;", strong("Data Statistics")),
            div(
              style = "text-align: justify;",
              "This page serves an ",
              strong("computation parameter overview."),
              br(),
              paste(
                "The statistics on the right will be generated considering",
                "the loaded dataset as well as the data settings provided."
              ),
              br(), br(),
              paste(
                "You will notice that changing certain settings in the",
                "Processed Data tab will change these statistics. This",
                "might be due to different settings in automatic annotation",
                "generation, where more or less datapoints will be ignored,",
                "or how certain data types are summarized."
              ),
              br(), br(),
              paste(
                "We recommend looking at these statistics carefully before",
                "performing Signature and Trajectory computation as part of",
                "the underlying problems that need to be solved for this",
                "computation are NP-complete and thus"
              ),
              strong("computationally complex."),
              paste(
                "The statistics given here represent parameters our",
                "computations scale with exponentially."
              )
            ),
            br(),
            div(
              style = "text-align: center;",
              actionButton(
                "generate_stats",
                "Generate Statistics",
                class = "btn btn-success",
                width = "100%"
              )
            )
          ),
          mainPanel(
            width = 9,
            div(
              style = "border-bottom: 1px solid lightgrey; padding: 10px;",
              h3("Cohort Statistics"),
              div(
                style = paste(
                  "display: flex; justify-content: space-between;",
                  "align-items: flex-start; margin-bottom: 10px;"
                ),
                div(
                  style = "width: 48%; text-align: left; font-weight: bold;",
                  "Number of Patients with X number of Clones"
                ),
                div(
                  style = "width: 48%; text-align: left; font-weight: bold;",
                  "Number of Patients with X number of Variants"
                )
              ),
              div(
                style = paste(
                  "display: flex; justify-content: space-between;",
                  "align-items: stretch;"
                ),
                div(
                  style = "width: 48%;",
                  plotOutput("clone_count_plot", height = "300px")
                ),
                div(
                  style = "width: 48%;",
                  plotOutput("variant_count_plot", height = "300px")
                )
              )
            ),
            div(
              style = "border-bottom: 1px solid lightgrey; padding: 10px;",
              h3("Variant Statistics"),
              div(
                style = paste(
                  "width: 100%; text-align: left; font-weight: bold;",
                  "margin-bottom: 10px;"
                ),
                paste(
                  "Distribution of Variants occuring in the different",
                  "chromosomes"
                )
              ),
              div(
                style = "width: 100%; padding: 10px; text-align: left;",
                plotOutput("chromosome_distribution_plot", height = "300px")
              ),
              div(
                style = paste(
                  "width: 100%; text-align: left; font-weight: bold;",
                  "margin-bottom: 10px;"
                ),
                "Distribution of Variants occuring in one/multiple patients"
              ),
              div(
                style = "width: 100%;",
                plotOutput("variants_multi_occurrence_plot", height = "300px")
              ),
              div(
                style = paste(
                  "display: flex; justify-content: space-between;",
                  "align-items: flex-start; margin-bottom: 10px;"
                ),
                div(
                  style = "width: 100%; text-align: left; font-weight: bold;",
                  "CNV to SNV ratio"
                )
              ),
              div(
                style = paste(
                  "display: flex; justify-content: space-between;",
                  "align-items: stretch; margin-bottom: 50px;"
                ),
                div(
                  style = "width: 100%;",
                  uiOutput("CNV_SNV_ratio_plot")
                )
              ),
              div(
                style = paste(
                  "display: flex; justify-content: space-between;",
                  "align-items: flex-start; margin-top: 50px;",
                  "margin-bottom: 10px;"
                ),
                div(
                  style = "width: 48%; text-align: left; font-weight: bold;",
                  "Most frequent Variants"
                )
              ),
              div(
                style = paste(
                  "display: flex; justify-content: space-between;",
                  "align-items: stretch; margin-bottom: 10px;"
                ),
                div(
                  style = "width: 100%;",
                  plotOutput("top_variants_plot", height = "500px")
                )
              )
            )
          )
        )
      ),
      tabPanel(
        "Signature Computation",
        sidebarLayout(
          sidebarPanel(
            width = 3,
            div(style = "font-size: 25px;", strong("Signature Computation")),
            div(
              style = "text-align: justify;",
              paste(
                "This tab is for the configuration of Signature computation",
                "and result showcasing. The settings below are for",
                "determining the computation methods and result in a",
                "complete dataframe including all signatures based on those",
                "settings. The settings in the top bar are applied",
                "afterwards in order to filter through the results",
                "efficiently. Changing the settings in the top bar will not",
                "require recomputation, changing the settings on the left",
                "will."
              ),
            ), tags$style(HTML("br { margin: 0; padding: 0; }")),
            div(
              style = "text-align: justify;",
              "Click on a computed signature if you are interested in more",
              "details.",
            ),
            br(),
            div(style = "font-size: 25px;", strong("Settings")),
            numericInput(
              "minimum_support_of_signatures",
              "Minimum number of occurences",
              value = 2,
              min = 1,
              max = 50
            ),
            numericInput(
              "signature_length_setter",
              "Max. Signature length",
              value = 2,
              min = 2,
              max = 5
            ),
            radioButtons(
              inputId = "radio_multiple_variants_per_clone",
              label = "Allow sequences within clone?",
              choices = c("Yes" = TRUE, "No" = FALSE),
              selected = TRUE,
              inline = TRUE
            ),
            # radioButtons(
            #   inputId = "ignore_single_occurrence",
            #   label = "Ignore single occurrence Variants? (TODO)",
            #   choices = c("Yes" = TRUE, "No" = FALSE),
            #   selected = TRUE,
            #   inline = TRUE
            # ),
            div(
              style = "text-align: center;",
              actionButton(
                "compute_signatures",
                "Compute Signatures",
                class = "btn btn-success",
                width = "100%"
              )
            ),
            br(),
            selectInput(
              inputId = "significance_estimation_method",
              label = "Method for p value estimation:",
              choices = c(
                "Uniform label probabilities" = "Sandmann",
                "Multiplicity-adjusted label probabilities" = "Sandmann Nick"
              ),
              selected = "Sandmann"
            ),
            div(
              style = "text-align: center;",
              actionButton(
                "estimate_significance",
                "Estimate Significance",
                class = "btn btn-success",
                width = "100%"
              )
            ),
            br(),
            strong("Optional CNV region analysis"),
            p("After computing signatures, compare shared regions of exact CNVs."),
            actionButton("postprocess_cnv_regions", "Analyze shared CNV regions",
              class = "btn btn-success", width = "100%"),
            div(
              style = "text-align: center; margin-top: 12px;",
              downloadButton("export_cnv_regions", "Export region signatures",
                class = "btn btn-success")
            )
          ),
          mainPanel(
            height = 10, width = 9,
            fluidRow(
              style = paste(
                "margin-bottom: 5px; background-color: #ecf0f1;",
                "padding-top: 5px; padding-right: 15px"
              ),
              column(
                3,
                textInput(
                  "signature_search",
                  "Search for specific variants:",
                  value = ""
                )
              ),
              column(
                2,
                numericInput(
                  "relative_frequency",
                  "Min. relative frequency:",
                  value = 0.000,
                  min = 0,
                  max = 1,
                  step = 0.001
                )
              ),
              column(
                2,
                numericInput(
                  "absolute_frequency_signatures",
                  "Min. occurences:",
                  value = 2,
                  min = 2,
                  step = 1
                )
              ),
              column(
                2,
                numericInput(
                  "min_p_value",
                  "Max. p value:",
                  value = 0.05,
                  min = 0,
                  max = 1,
                  step = 0.001
                )
              ),
              column(
                3,
                div(
                  style = "text-align: center; padding-top: 25px;",
                  actionButton(
                    "apply_signature_filter",
                    "Enable Filter",
                    class = "btn btn-success",
                    width = "100%"
                  )
                )
              )
            ),
            style = "padding: 0;",
            uiOutput("signature_computation_main_panel"),
            uiOutput("cnv_region_panel")
          )
        )
      ),
      tabPanel(
        "Trajectory Approximation",
        sidebarLayout(
          sidebarPanel(
            width = 3,
            div(
              style = "text-align: justify;",
              div(
                style = "font-size: 25px;",
                strong("Trajectory approximation")
              ),
              paste(
                "Previously computed Signatures (linear patterns) can be",
                "combined to identify significantly enriched branching",
                "structures, which we refer to as trajectories. By pressing",
                "the"
              ),
              em("Compute trajectories"),
              paste(
                "button, the program will automatically filter all",
                "Signatures given the filter configuration in the"
              ),
              em("Signature Computation"),
              paste(
                "tab, no matter if the filter view is enabled or disabled."
              )
            ), br(),
            numericInput(
              "minimum_support_of_trajectories",
              "Minimum number of occurences",
              value = 2,
              min = 2,
              max = 50
            ),
            numericInput(
              "recombination_number_setter",
              "Number of recombinations",
              value = 1,
              min = 1,
              max = 10
            ),
            radioButtons(
              inputId = "compute_independent_root_branches",
              label = "Include independent branches through the root?",
              choices = c("Yes" = TRUE, "No" = FALSE),
              selected = TRUE,
              inline = TRUE
            ),
            radioButtons(
              inputId = "multiple_identical_variants_for_branches",
              label = "Allow the same variant to form multiple subclones?",
              choices = c("Yes" = TRUE, "No" = FALSE),
              selected = FALSE,
              inline = TRUE
            ),
            div(
              style = "text-align: center;",
              actionButton(
                "compute_trajectories",
                "Compute trajectories",
                class = "btn btn-success",
                width = "100%"
              )
            ), br(),
            selectInput(
              inputId = "significance_estimation_method_trajectories",
              label = "Method for p value estimation:",
              choices = c("Uniform label probabilities" = "Sandmann"),
              selected = "Sandmann"
            ),
            selectInput(
              inputId = "p_value_adjustment_trajectories",
              label = "p value adjustment:",
              choices = c("None", "Benjamini-Hochberg"),
              selected = "None"
            ),
            # selectInput(
            #   inputId = "aposteriori_cnv_refactoring",
            #   label = paste(
            #     "Refactor CNVs?",
            #     "(Only applicable for exact CNV encoding)"
            #   ),
            #   choices = c("No", "Interpret deletions as Gene mutations"),
            #   selected = "No"
            # ),
            div(
              style = "text-align: center;",
              actionButton(
                "estimate_significant_trajectories",
                "Estimate significant trajectories",
                class = "btn btn-success",
                width = "100%"
              )
            ), br(),
            div(
              style = "text-align: center;",
              actionButton(
                "export_trajectories",
                "Export Trajectory Data",
                class = "btn btn-success",
                width = "100%"
              )
            )
          ),
          mainPanel(
            height = 10, width = 9,
            fluidRow(
              style = paste(
                "margin-bottom: 5px; background-color: #ecf0f1;",
                "padding-top: 5px; padding-right: 15px"
              ),
              column(
                3,
                textInput(
                  "variant_search_trajectories",
                  "Search for specific variants:",
                  value = ""
                )
              ),
              column(
                2,
                numericInput(
                  "relative_frequency_trajectories",
                  "Min. relative frequency:",
                  value = 0.000,
                  min = 0,
                  max = 1,
                  step = 0.001
                )
              ),
              column(
                2,
                numericInput(
                  "absolute_frequency_trajectories",
                  "Min. occurences:",
                  value = 2,
                  min = 2,
                  step = 1
                )
              ),
              column(
                2,
                numericInput(
                  "min_p_value_trajectories",
                  "Max. p value:",
                  value = 1.000,
                  min = 0,
                  max = 1,
                  step = 0.001
                )
              ),
              column(
                3,
                div(
                  style = "text-align: center; padding-top: 25px;",
                  actionButton(
                    "apply_trajectory_filter",
                    "Enable Filter",
                    class = "btn btn-success",
                    width = "100%"
                  )
                )
              )
            ),
            style = "padding: 0;",
            uiOutput("trajectory_approximation_main_panel")
          )
        )
      )
    ),
    tags$head(
      tags$style(HTML("
      .sidebar-btn {
        display: block;
        width: 100%;
        margin-bottom: 10px;
        height: 40px;
        text-align: center;
        font-size: 14px;
        padding: 5px;
        box-sizing: border-box;
      }

      .sidebar-btn .glyphicon {
        display: none;
      }
    "))
    )
  )
}
