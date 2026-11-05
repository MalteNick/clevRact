handle_loaded_dataframe <- function(state, output, input) {
  tryCatch(
    {
      state$data_set <- read.csv2(
        input$file$datapath,
        fileEncoding = "ISO-8859-1",
        sep = ";"
      )
      required_columns <- c(
        "Patient.ID", "Variant", "Clone", "Parent", "Variant.Type",
        "Chromosome", "Start", "End"
      )
      missing_columns <- setdiff(required_columns, colnames(state$data_set))
      if (length(missing_columns) > 0) {
        shinyalert(
          text = paste(
            "The following required columns are missing:",
            paste(missing_columns, collapse = ", ")
          ),
          type = "error",
          showCancelButton = FALSE,
          confirmButtonText = "OK",
          closeOnEsc = TRUE,
          closeOnClickOutside = FALSE
        )
        return()
      }
      ti_columns <- grep("^T[0-9]+$", colnames(state$data_set), value = TRUE)
      if (length(ti_columns) > 0) {
        ti_numbers <- as.integer(sub("T", "", ti_columns))
        missing_ti <- setdiff(seq_len(max(ti_numbers)), ti_numbers)
        if (length(missing_ti) > 0) {
          shinyalert(
            title = "Missing T Columns",
            text = paste("Missing T columns: T", min(missing_ti)),
            type = "error",
            showCancelButton = FALSE,
            confirmButtonText = "OK",
            closeOnEsc = TRUE,
            closeOnClickOutside = FALSE
          )
          return()
        }
      }
      # Value correction
      state$data_set <- state$data_set[
        ,
        union(required_columns, ti_columns),
        drop = FALSE
      ]
      state$data_set$Variant.Type <- ifelse(
        state$data_set$Variant.Type %in% c(
          "CNV_del", "CNV_dup", "CNV_loh", "SNV/INDEL", "NA"
        ),
        state$data_set$Variant.Type,
        "NA"
      )
      for (col in ti_columns) {
        state$data_set[[col]] <- ifelse(
          state$data_set[[col]] %in% c(TRUE, FALSE),
          state$data_set[[col]],
          TRUE
        )
      }
      state$patient_ids <- unique(state$data_set$Patient.ID)
      state$current_patient_id <- state$patient_ids[[1]]
      state$csv_loaded <- TRUE
      output$id_selector <- renderUI({
        selectizeInput("id_select",
          "Select Patient ID",
          choices = state$patient_ids,
          selected = state$current_patient_id
        )
      })
      state$dataframe_edited <- FALSE
      generate_editor_main_panel(state, input, output)
      return(TRUE)
    },
    error = function(e) {
      shinyalert(
        title = "Error Reading Dataset",
        text = paste("There was an error reading the dataset:", e$message),
        type = "error",
        showCancelButton = FALSE,
        confirmButtonText = "OK",
        closeOnEsc = TRUE,
        closeOnClickOutside = FALSE
      )
      return(FALSE)
    }
  )
}

generate_editor_main_panel <- function(state, input, output) {
  if (
    is.null(state$data_set) ||
      !state$csv_loaded ||
      is.null(state$current_patient_id)
  ) {
    output$input_data_main_panel <- renderUI({
      div(
        img(
          src = "clevRact-assets/clevRact_bw.png",
          height = "150px",
          style = "margin-bottom: 10px;"
        ),
        div(
          HTML(
            paste0(
              "No dataset generated. Load a dataset or create one by ",
              "clicking the&nbsp;",
              as.character(em("New Dataset")),
              "&nbsp;button."
            )
          ),
          style = "white-space: nowrap; max-width: 80%;"
        ),
        style = "
          display: flex;
          flex-direction: column;
          justify-content: center;
          align-items: center;
          height: 60vh;
          text-align: center;
        "
      )
    })
    output$data_table <- NULL
    return()
  }
  if (nrow(state$data_set) == 0) {
    output$input_data_main_panel <- renderUI({
      div(
        "Empty dataset. Add a patient to populate it.",
        style = "
          display: flex;
          justify-content: center;
          align-items: center;
          height: 60vh;
          text-align: center;
        "
      )
    })
    output$data_table <- NULL
    return()
  }
  if (!input$toggle_view) {
    output$input_data_main_panel <- renderUI({
      DT::dataTableOutput("data_table")
    })
    output$data_table <- DT::renderDataTable({
      current_data <- state$data_set[
        state$data_set$Patient.ID == state$current_patient_id, ,
        drop = FALSE
      ]
      if (state$dataset_selection == "none") {
        editable_boolean <- TRUE
      } else {
        editable_boolean <- FALSE
      }
      current_data <- current_data[
        ,
        c(
          c(
            "Patient.ID", "Variant", "Clone", "Parent", "Chromosome",
            "Start", "End", "Variant.Type"
          ),
          grep("^T", colnames(state$data_set), value = TRUE)
        ),
        drop = FALSE
      ]
      datatable(current_data,
        options = list(
          scrollX = TRUE,
          dom = "t",
          paging = FALSE,
          searching = FALSE,
          scrollCollapse = TRUE,
          fixedHeader = TRUE,
          columnDefs = list(
            list(className = "dt-left", targets = "_all"),
            list(visible = FALSE, targets = "Patient.ID")
          ),
          order = list(list(
            which(colnames(current_data) == "Clone") - 1,
            "asc"
          ))
        ),
        editable = editable_boolean,
        rownames = FALSE,
        colnames = c(
          "Variant Type" = "Variant.Type",
          "Patient ID" = "Patient.ID"
        ),
        selection = state$dataset_selection
      )
    })
  } else {
    seaObject <- NULL
    errorOccurred <- FALSE
    tryCatch(
      {
        seaObject <- createSeaObjectFromDataset(
          state$data_set[
            state$data_set$Patient.ID == state$current_patient_id, ,
            drop = FALSE
          ]
        )
      },
      error = function(e) {
        errorOccurred <- TRUE
        output$input_data_main_panel <- renderUI({
          div(
            paste("Provided patient data cannot be processed:", e$message),
            style = paste(
              "display: flex; justify-content: center;",
              "align-items: center; height: 60vh;",
              "text-align: center;"
            )
          )
        })
      }
    )
    if (errorOccurred || is.null(seaObject)) {
      return()
    }
    output$shark_plot <- renderPlot({
      sharkPlot(seaObject, showLegend = TRUE)
    })
    output$dolphin_plot <- renderPlot({
      dolphinPlot(seaObject, showLegend = TRUE)
    })
    output$input_data_main_panel <- renderUI({
      column(
        width = 12,
        fluidRow(
          column(
            width = 6,
            div(
              style = "text-align: center; font-weight: bold;",
              "Dolphin Plot"
            )
          ),
          column(
            width = 6,
            div(style = "text-align: center; font-weight: bold;", "Shark Plot")
          )
        ),
        fluidRow(
          column(
            width = 6,
            plotOutput("dolphin_plot")
          ),
          column(
            width = 6,
            plotOutput("shark_plot")
          )
        )
      )
    })
  }
}

generate_processed_data_main_panel <- function(state, output) {
  if (
    is.null(state$preprocessed_data) ||
      is.null(state$current_processed_patient_id)
  ) {
    output$processed_data_main_panel <- renderUI({
      div(
        img(
          src = "clevRact-assets/clevRact_bw.png",
          height = "150px",
          style = "margin-bottom: 10px;"
        ),
        div(
          HTML(
            paste0(
              "No dataset processed. Configure your settings and click the ",
              em("Generate Processed Data"),
              " button."
            )
          ),
          style = "white-space: nowrap; max-width: 80%;"
        ),
        style = "
          display: flex;
          flex-direction: column;
          justify-content: center;
          align-items: center;
          height: 60vh;
          text-align: center;
        "
      )
    })
    return()
  }
  output$processed_data_main_panel <- renderUI({
    DT::dataTableOutput("processed_data_table")
  })
  output$processed_data_table <- DT::renderDataTable({
    current_data <- state$preprocessed_data[[1]][
      state$preprocessed_data[[1]]$Patient.ID ==
        state$current_processed_patient_id, ,
      drop = FALSE
    ]
    datatable(current_data,
      options = list(
        scrollX = TRUE,
        dom = "t",
        paging = FALSE,
        searching = FALSE,
        scrollCollapse = TRUE,
        fixedHeader = TRUE,
        columnDefs = list(
          list(className = "dt-left", targets = "_all"),
          list(visible = FALSE, targets = c("Patient.ID", "Chromosome"))
        )
      ),
      editable = FALSE,
      rownames = FALSE,
      colnames = c(
        "Variant Type" = "Variant.Type",
        "Patient ID" = "Patient.ID"
      ),
      selection = "none"
    )
  })
}

##########################################
# Functions to reset tabs to default states
##########################################

set_input_data_panel_to_blank <- function(output, state) {
  output$input_data_main_panel <- renderUI({
    div(
      img(
        src = "clevRact-assets/clevRact_bw.png",
        height = "150px",
        style = "margin-bottom: 10px;"
      ),
      div(
        HTML(
            paste0(
            "No dataset generated. Load a dataset or create one by ",
            "clicking the&nbsp;",
            as.character(em("New Dataset")),
            "&nbsp;button."
          )
        ),
        style = "white-space: nowrap; max-width: 80%;"
      ),
      style = "
        display: flex;
        flex-direction: column;
        justify-content: center;
        align-items: center;
        height: 60vh;
        text-align: center;
      "
    )
  })
  output$id_selector <- renderUI({
    selectizeInput("id_select",
      "Select Patient ID",
      choices = list()
    )
  })
}

set_processed_data_panel_to_blank <- function(output, state) {
  output$processed_id_selector <- renderUI({
    selectizeInput("processed_id_select",
      "Select Patient ID",
      choices = list()
    )
  })
  output$processed_data_main_panel <- renderUI({
    div(
      img(
        src = "clevRact-assets/clevRact_bw.png",
        height = "150px",
        style = "margin-bottom: 10px;"
      ),
      div(
        HTML(
          paste0(
            "No dataset processed. Configure your settings and click the ",
            as.character(em("Generate Processed Data")),
            "&nbsp;button."
          )
        ),
        style = "white-space: nowrap; max-width: 80%;"
      ),
      style = "
        display: flex;
        flex-direction: column;
        justify-content: center;
        align-items: center;
        height: 60vh;
        text-align: center;
      "
    )
  })
  output$processed_data_table <- NULL
  state$current_processed_patient_id <- NULL
  state$preprocessed_data <- NULL
}

set_data_statistics_panel_to_blank <- function(output, state) {
  output$clone_count_plot <- renderPlot({
    ggplot() +
      annotate(
        "text",
        x = 0,
        y = 10,
        label = "No statistics generated",
        size = 6
      ) +
      scale_x_continuous(expand = c(0, 0), breaks = seq(0, 10, by = 1)) +
      scale_y_continuous(expand = c(0, 0), breaks = seq(0, 10, by = 1)) +
      theme_bw() +
      theme(
        axis.text = element_text(size = 12),
        axis.title = element_blank()
      )
  })

  output$variant_count_plot <- renderPlot({
    ggplot() +
      annotate(
        "text",
        x = 0,
        y = 10,
        label = "No statistics generated",
        size = 6
      ) +
      scale_x_continuous(expand = c(0, 0), breaks = seq(0, 10, by = 1)) +
      scale_y_continuous(expand = c(0, 0), breaks = seq(0, 10, by = 1)) +
      theme_bw() +
      theme(
        axis.text = element_text(size = 12),
        axis.title = element_blank()
      )
  })

  output$variants_multi_occurrence_plot <- renderPlot({
    ggplot() +
      annotate(
        "text",
        x = 0,
        y = 10,
        label = "No statistics generated",
        size = 6
      ) +
      scale_x_continuous(expand = c(0, 0), breaks = seq(0, 10, by = 1)) +
      scale_y_continuous(expand = c(0, 0), breaks = seq(0, 10, by = 1)) +
      theme_bw() +
      theme(
        axis.text = element_text(size = 12),
        axis.title = element_blank()
      )
  })

  output$top_variants_plot <- renderPlot({
    ggplot() +
      annotate(
        "text",
        x = 0,
        y = 10,
        label = "No statistics generated",
        size = 6
      ) +
      scale_x_continuous(expand = c(0, 0), breaks = seq(0, 10, by = 1)) +
      scale_y_continuous(expand = c(0, 0), breaks = seq(0, 10, by = 1)) +
      theme_bw() +
      theme(
        axis.text = element_text(size = 12),
        axis.title = element_blank()
      )
  })

  output$CNV_SNV_ratio_plot <- renderUI({
    percentage <- 30
    snv_width <- paste0(percentage, "%")
    cnv_width <- paste0(100 - percentage, "%")
    div(
      style = paste(
        "width: 100%; height: 30px;",
        "display: flex; border: 1px solid black;"
      ),
      div(
        style = paste0(
          "width:", snv_width,
          "; height: 100%; background-color: white;"
        )
      )
    )
  })

  output$chromosome_distribution_plot <- renderPlot({
    ggplot() +
      annotate(
        "text",
        x = 0,
        y = 10,
        label = "No statistics generated",
        size = 6
      ) +
      scale_x_continuous(expand = c(0, 0), breaks = seq(0, 10, by = 1)) +
      scale_y_continuous(expand = c(0, 0), breaks = seq(0, 10, by = 1)) +
      theme_bw() +
      theme(
        axis.text = element_text(size = 12),
        axis.title = element_blank()
      )
  })
}

set_signature_computation_panel_to_blank <- function(output, state) {
  output$signature_computation_main_panel <- renderUI({
    div(
      img(
        src = "clevRact-assets/clevRact_bw.png",
        height = "150px",
        style = "margin-bottom: 10px;"
      ),
      div(
        HTML(
          paste0(
            "No Signatures computed. Configure your settings and click the ",
            as.character(em("Compute Signatures")),
            "&nbsp;button."
          )
        ),
        style = "white-space: nowrap; max-width: 80%;"
      ),
      style = "
        display: flex;
        flex-direction: column;
        justify-content: center;
        align-items: center;
        height: 60vh;
        text-align: center;
      "
    )
  })
  output$signature_table <- NULL
  state$signature_data_set <- NULL
  state$filtered_signature_data_set <- NULL
  state$signature_side_information <- NULL
  state$signature_type <- 0
  state$cnv_region_result <- NULL
  output$cnv_region_panel <- NULL
}

set_trajectory_approximation_panel_to_blank <- function(output, state) {
  output$trajectory_approximation_main_panel <- renderUI({
    div(
      img(
        src = "clevRact-assets/clevRact_bw.png",
        height = "150px",
        style = "margin-bottom: 10px;"
      ),
      div(
        HTML(
          paste0(
            "No trajectories computed. Configure your settings and click the ",
            as.character(em("Compute Trajectories")),
            "&nbsp;button."
          )
        ),
        style = "white-space: nowrap; max-width: 80%;"
      ),
      style = "
        display: flex;
        flex-direction: column;
        justify-content: center;
        align-items: center;
        height: 60vh;
        text-align: center;
      "
    )
  })
  output$trajectory_table <- NULL
  state$trajectory_data_set <- NULL
  state$filtered_trajectory_data_set <- NULL
}

.build_server <- function() {
  function(input, output, session) {
    state <- reactiveValues(
      csv_loaded = FALSE,
      dataframe_edited = FALSE,
      data_set = NULL,
      patient_ids = NULL,
      current_patient_id = NULL,
      preprocessed_data = NULL,
      current_processed_patient_id = NULL,
      signature_data_set = NULL,
      filtered_signature_data_set = NULL,
      signature_side_information = NULL,
      signature_type = 0,
      trajectory_data_set = NULL,
      filtered_trajectory_data_set = NULL,
      variant_edit_mode = "enable selection",
      dataset_selection = "none"
    )
    set_input_data_panel_to_blank(output, state)
    set_processed_data_panel_to_blank(output, state)
    set_data_statistics_panel_to_blank(output, state)
    set_signature_computation_panel_to_blank(output, state)
    set_trajectory_approximation_panel_to_blank(output, state)

    ######################################
    # Eventhandlers in the data editor tab
    ######################################

    observeEvent(input$new_dataset, {
      if (!is.null(state$data_set) && state$dataframe_edited) {
        shinyalert(
          title = "Unsaved Changes",
          text = paste(
            "There have been edits to this dataset that cannot be",
            "tracked in the original file. To save the changes click the",
            "Save Dataset button first to download the dataset including",
            "all changes. Are you sure you want to load a new dataset?"
          ),
          type = "warning",
          showCancelButton = TRUE,
          confirmButtonText = "Discard Changes & Load",
          cancelButtonText = "Cancel",
          closeOnEsc = TRUE,
          closeOnClickOutside = FALSE,
          callbackR = function(confirm) {
            if (confirm) {
              showModal(modalDialog(
                title = "No Dataset Loaded",
                numericInput(
                  "num_timepoints",
                  "Number of timepoints (0 for none):",
                  min = 0,
                  value = 0
                ),
                footer = tagList(
                  modalButton("Cancel"),
                  actionButton("confirm_create_dataset", "Create New Dataset")
                ),
                easyClose = TRUE
              ))
            }
          }
        )
      } else {
        showModal(modalDialog(
          title = "No Dataset Loaded",
          numericInput(
            "num_timepoints",
            "Number of timepoints (0 for none):",
            min = 0,
            value = 0
          ),
          footer = tagList(
            modalButton("Cancel"),
            actionButton("confirm_create_dataset", "Create New Dataset")
          ),
          easyClose = TRUE
        ))
      }
    })

    observeEvent(input$load_dataset, {
      if (!is.null(state$data_set)) {
        if (state$dataframe_edited) {
          shinyalert(
            title = "Unsaved Changes",
            text = paste(
              "There have been edits to this dataset that cannot be",
              "tracked in the original file. To save the changes click the",
              "Save Dataset button first to download the dataset including",
              "all changes. Are you sure you want to load a new dataset?"
            ),
            type = "warning",
            showCancelButton = TRUE,
            confirmButtonText = "Discard Changes & Load",
            cancelButtonText = "Cancel",
            closeOnEsc = TRUE,
            closeOnClickOutside = FALSE,
            callbackR = function(confirm) {
              if (confirm) {
                showModal(modalDialog(
                  title = "Select a file",
                  fileInput("file", "Choose a file", accept = c(".csv")),
                  easyClose = TRUE
                ))
              }
            }
          )
        } else {
          showModal(modalDialog(
            title = "Select a file",
            fileInput("file", "Choose a file", accept = c(".csv")),
            easyClose = TRUE
          ))
        }
      } else {
        showModal(modalDialog(
          title = "Select a file",
          fileInput("file", "Choose a file", accept = c(".csv")),
          easyClose = TRUE
        ))
      }
    })

    observeEvent(input$id_select, {
      state$current_patient_id <- input$id_select
    })

    observeEvent(input$next_patient, {
      if (!is.null(state$data_set)) {
        current_index <- which(state$patient_ids == state$current_patient_id)
        next_index <- ifelse(current_index == length(state$patient_ids),
          1, current_index + 1
        )
        state$current_patient_id <- state$patient_ids[next_index]
        output$id_selector <- renderUI({
          selectizeInput("id_select",
            "Select Patient ID",
            choices = state$patient_ids,
            selected = state$current_patient_id
          )
        })
      } else {
        showNotification(
          "Cannot change patient without a dataset",
          type = "message",
          duration = 3
        )
      }
    })

    observeEvent(input$previous_patient, {
      if (!is.null(state$data_set)) {
        current_index <- which(state$patient_ids == state$current_patient_id)
        prev_index <- ifelse(
          current_index == 1,
          length(state$patient_ids),
          current_index - 1
        )
        state$current_patient_id <- state$patient_ids[prev_index]
        output$id_selector <- renderUI({
          selectizeInput("id_select",
            "Select Patient ID",
            choices = state$patient_ids,
            selected = state$current_patient_id
          )
        })
      } else {
        showNotification(
          "Cannot change patient without a dataset",
          type = "message",
          duration = 3
        )
      }
    })

    observeEvent(input$change_timepoints, {
      if (!is.null(state$data_set)) {
        showModal(modalDialog(
          title = "Change Timepoints",
          numericInput(
            "time_point_input",
            "Please enter the number of timepoints for the dataset:",
            value = length(grep(
              "^T",
              colnames(state$data_set),
              value = TRUE
            )),
            min = 0,
            step = 1
          ),
          footer = tagList(
            modalButton("Cancel"),
            actionButton("confirm_timepoint_change", "Confirm")
          ),
          size = "m"
        ))
      } else {
        showNotification(
          "Cannot change timepoints without a dataset",
          type = "message",
          duration = 3
        )
      }
    })

    observeEvent(input$confirm_timepoint_change, {
      k <- length(grep("^T", colnames(state$data_set), value = TRUE))
      if (input$time_point_input > k) {
        num_columns_to_add <- input$time_point_input - k
        new_columns <- paste0("T", seq.int(k + 1, k + num_columns_to_add))
        state$data_set[new_columns] <- "TRUE"
        showNotification(
          paste(
            "Added",
            num_columns_to_add,
            "new Timepoint columns.",
            collapse = " "
          ),
          type = "message",
          duration = 3
        )
      } else if (input$time_point_input < k) {
        num_columns_to_remove <- k - input$time_point_input
        columns_to_remove <- tail(
          grep("^T", colnames(state$data_set), value = TRUE),
          num_columns_to_remove
        )
        state$data_set <- state$data_set[
          ,
          !colnames(state$data_set) %in% columns_to_remove
        ]
        showNotification(
          paste(
            "Removed",
            num_columns_to_remove,
            "Timepoint columns.",
            collapse = " "
          ),
          type = "message",
          duration = 3
        )
      }
      generate_editor_main_panel(state, input, output)
      removeModal()
    })

    observeEvent(input$add_patient, {
      if (!is.null(state$data_set)) {
        showModal(modalDialog(
          title = "Add New Patient",
          textInput("new_patient_id", "Enter new patient ID:", ""),
          footer = tagList(
            modalButton("Cancel"),
            actionButton("confirm_add_patient", "Add")
          ),
          easyClose = TRUE
        ))
        state$dataframe_edited <- TRUE
      } else {
        showNotification(
          "Cannot add patient without a dataset",
          type = "message",
          duration = 3
        )
      }
    })

    observeEvent(input$confirm_create_dataset, {
      num_timepoints <- input$num_timepoints
      if (is.null(num_timepoints) || num_timepoints < 0) {
        showModal(modalDialog(
          title = "Invalid Input",
          "Please specify a valid number of timepoints (0 or more).",
          footer = modalButton("OK"),
          easyClose = TRUE
        ))
        return()
      }

      required_columns <- c(
        "Patient.ID", "Variant", "Clone", "Parent", "Variant.Type",
        "Chromosome", "Start", "End"
      )
      ti_columns <- if (num_timepoints > 0) {
        paste("T", seq_len(num_timepoints), sep = "")
      } else {
        character(0)
      }
      all_columns <- c(required_columns, ti_columns)

      state$data_set <- as.data.frame(matrix(
        nrow = 0,
        ncol = length(all_columns),
        dimnames = list(NULL, all_columns)
      ), stringsAsFactors = FALSE)

      state$patient_ids <- unique(state$data_set$Patient.ID)
      state$csv_loaded <- TRUE

      output$id_selector <- renderUI({
        selectizeInput("id_select",
          "Select Patient ID",
          choices = list()
        )
      })
      generate_editor_main_panel(state, input, output)
      state$dataframe_edited <- TRUE
      removeModal()
      showModal(modalDialog(
        title = "New Dataset Created",
        "You can now add a new patient.",
        footer = tagList(
          modalButton("OK"),
          actionButton("add_patient_after_create", "Add Patient")
        ),
        easyClose = TRUE
      ))
    })

    observeEvent(input$add_patient_after_create, {
      removeModal()
      showModal(modalDialog(
        title = "Add New Patient",
        textInput("new_patient_id", "Enter new patient ID:", ""),
        footer = tagList(
          modalButton("Cancel"),
          actionButton("confirm_add_patient", "Add")
        ),
        easyClose = TRUE
      ))
    })

    observeEvent(input$confirm_add_patient,
      {
        new_id <- input$new_patient_id
        if (is.null(new_id) || new_id == "" || new_id %in% state$patient_ids) {
          showModal(modalDialog(
            title = "Error",
            if (is.null(new_id) || new_id == "") {
              "Patient ID cannot be empty."
            } else {
              paste("Patient ID", new_id, "already exists.")
            },
            footer = modalButton("OK"),
            easyClose = TRUE
          ))
          return()
        }
        new_row <- as.data.frame(matrix(
          nrow = 1,
          ncol = ncol(state$data_set),
          dimnames = list(NULL, colnames(state$data_set))
        ), stringsAsFactors = FALSE)
        new_row[["Patient.ID"]] <- new_id
        state$data_set <- rbind(state$data_set, new_row)
        state$patient_ids <- unique(state$data_set$Patient.ID)
        output$id_selector <- renderUI({
          selectizeInput("id_select",
            "Select Patient ID",
            choices = state$patient_ids,
            selected = new_id
          )
        })
        state$current_patient_id <- new_id
        removeModal()
        state$dataframe_edited <- TRUE
        generate_editor_main_panel(state, input, output)
      },
      ignoreInit = TRUE
    )

    observeEvent(input$delete_patient, {
      if (is.null(state$data_set) || !state$csv_loaded) {
        showNotification(
          "Cannot delete patient without a dataset",
          type = "message",
          duration = 3
        )
        return()
      }
      if (is.null(state$current_patient_id)) {
        return()
      }
      shinyalert(
        title = "Confirm Deletion",
        text = paste(
          "Are you sure you want to delete this Patient ID",
          state$current_patient_id,
          "?"
        ),
        type = "warning",
        showCancelButton = TRUE,
        confirmButtonText = "Delete",
        cancelButtonText = "Cancel",
        closeOnEsc = TRUE,
        closeOnClickOutside = FALSE,
        callbackR = function(confirm) {
          if (confirm) {
            state$data_set <- state$data_set[
              state$data_set$Patient.ID != state$current_patient_id, ,
              drop = FALSE
            ]
            state$patient_ids <- unique(state$data_set$Patient.ID)
            if (length(state$patient_ids) > 0) {
              current_index <- which(
                state$patient_ids == state$current_patient_id
              )
              if (length(current_index) == 0 || current_index == 1) {
                state$current_patient_id <- state$patient_ids[1]
              } else {
                state$current_patient_id <- state$patient_ids[
                  current_index - 1
                ]
              }
            } else {
              state$current_patient_id <- NULL
            }
            state$dataframe_edited <- TRUE
            generate_editor_main_panel(state, input, output)
          }
        }
      )
    })

    observeEvent(input$add_variant, {
      if (is.null(state$data_set) || !state$csv_loaded) {
        showNotification(
          "Cannot add variant without a dataset",
          type = "message",
          duration = 3
        )
        return()
      }
      dataset_columns <- colnames(state$data_set)
      new_variant <- data.frame(
        Patient.ID = state$current_patient_id,
        stringsAsFactors = FALSE
      )
      for (col in dataset_columns) {
        if (col != "Patient.ID") {
          new_variant[[col]] <- NA
        }
      }
      state$data_set <- rbind(state$data_set, new_variant)
      state$dataframe_edited <- TRUE
      generate_editor_main_panel(state, input, output)
    })

    observeEvent(input$enable_disable_selection, {
      if (state$variant_edit_mode == "enable selection") {
        updateActionButton(
          session,
          "enable_disable_selection",
          label = "Untoggle Variant Selection"
        )
        state$variant_edit_mode <- "disable selection"
        state$dataset_selection <- "multiple"
        generate_editor_main_panel(state, input, output)
        # TODO bottom right popup to tell the user to load dataset
      } else if (state$variant_edit_mode == "disable selection") {
        updateActionButton(
          session,
          "enable_disable_selection",
          label = "Toggle Variant Selection"
        )
        state$variant_edit_mode <- "enable selection"
        state$dataset_selection <- "none"
        generate_editor_main_panel(state, input, output)
      }
    })

    observeEvent(input$delete_selection, {
      if (is.null(state$data_set)) {
        showNotification(
          "Cannot delete variants without a dataset",
          type = "message",
          duration = 3
        )
        return()
      }
      selected_rows <- input$data_table_rows_selected
      if (length(selected_rows) > 0) {
        patient_data <- state$data_set[
          state$data_set$Patient.ID == state$current_patient_id, ,
          drop = FALSE
        ]
        global_rows_to_delete <- which(
          state$data_set$Patient.ID == state$current_patient_id
        )[selected_rows]
        state$data_set <- state$data_set[-global_rows_to_delete, ]
        state$dataframe_edited <- TRUE
        generate_editor_main_panel(state, input, output)
      }
    })

    observeEvent(input$toggle_view, {
      generate_editor_main_panel(state, input, output)
    })

    observeEvent(input$data_table_cell_edit, {
      info <- input$data_table_cell_edit
      row <- info$row
      col <- info$col + 1
      new_value <- info$value
      fixed_columns <- c(
        "Patient.ID", "Variant", "Clone", "Parent", "Chromosome", "Start",
        "End", "Variant.Type"
      )
      col_name <- c(
        fixed_columns,
        grep("^T", colnames(state$data_set), value = TRUE)
      )[col]
      numeric_columns <- c("Clone", "Parent", value = TRUE)
      # Validate Variant Type Input
      if (col_name == "Variant.Type") {
        allowed_variant_types <- c(
          "CNV_del", "CNV_dup", "CNV_loh", "SNV/INDEL", "NA"
        )
        if (!(new_value %in% allowed_variant_types)) {
          showNotification(
            paste(
              "Invalid input. Please enter either 'CNV_del', 'CNV_dup',",
              "'CNV_loh', 'SNV/INDEL' or 'NA'"
            ),
            type = "error",
            duration = 5
          )
          generate_editor_main_panel(state, input, output)
          return()
        }
      }
      if (col_name == "Chromosome") {
        if (new_value == "x") {
          new_value <- "X"
        } else if (new_value == "y") {
          new_value <- "Y"
        }
        allowed_chromosomes <- c(as.character(seq_len(22)), "X", "Y")
        if (!(as.character(new_value) %in% allowed_chromosomes)) {
          showNotification("Invalid input. Please enter either 1-22 or X/Y",
            type = "error",
            duration = 5
          )
          generate_editor_main_panel(state, input, output)
          return()
        }
      }
      if (col_name == "Clone" || col_name == "Parent") {
        if (new_value < 0) {
          showNotification("Invalid input. Please enter a non-negative value.",
            type = "error",
            duration = 5
          )
          generate_editor_main_panel(state, input, output)
          return()
        }
      }
      if (col_name %in% numeric_columns) {
        if (is.na(as.numeric(new_value))) {
          showNotification(
            "Invalid input. Please enter a numeric value.",
            type = "error"
          )
          generate_editor_main_panel(state, input, output)
          return()
        }
      }
      if (grepl("^T\\d+$", col_name)) {
        if (!(new_value %in% c("TRUE", "FALSE", TRUE, FALSE))) {
          showNotification(
            "Invalid input. Please enter either TRUE or FALSE",
            type = "error",
            duration = 5
          )
          generate_editor_main_panel(state, input, output)
          return()
        }

        logical_value <- as.logical(new_value)

        patient_data <- state$data_set[
          state$data_set$Patient.ID == state$current_patient_id, ,
          drop = FALSE
        ]
        row_index <- which(
          state$data_set$Patient.ID == state$current_patient_id
        )[row]
        clone_value <- state$data_set[row_index, "Clone"]

        matching_rows <- which(
          state$data_set$Patient.ID == state$current_patient_id &
            state$data_set$Clone == clone_value
        )

        state$data_set[matching_rows, col_name] <- logical_value
      } else if (grepl("CCF", col_name, fixed = TRUE)) {
        patient_data <- state$data_set[
          state$data_set$Patient.ID == state$current_patient_id, ,
          drop = FALSE
        ]
        row_index <- which(
          state$data_set$Patient.ID == state$current_patient_id
        )[row]
        clone_value <- state$data_set[row_index, "Clone"]
        matching_rows <- which(
          state$data_set$Patient.ID == state$current_patient_id &
            state$data_set$Clone == clone_value
        )
        state$data_set[matching_rows, col_name] <- as.numeric(new_value)
      } else {
        patient_data <- state$data_set[
          state$data_set$Patient.ID == state$current_patient_id, ,
          drop = FALSE
        ]
        row_index <- which(
          state$data_set$Patient.ID == state$current_patient_id
        )[row]
        state$data_set[row_index, col_name] <- new_value
      }

      state$dataframe_edited <- TRUE
    })

    ######################################
    # Eventhandlers in the processed data tab
    ######################################

    observeEvent(input$processed_id_select, {
      state$current_processed_patient_id <- input$processed_id_select
    })

    observeEvent(input$dropdown_variant_generation, {
      if (input$dropdown_variant_generation == "No") {
        disable("dropdown_summary_snvs")
        disable("dropdown_summary_cnvs")
      } else {
        enable("dropdown_summary_snvs")
        enable("dropdown_summary_cnvs")
      }
    })

    observeEvent(input$next_processed_patient, {
      processed_patient_ids <- unique(state$preprocessed_data[[1]]$Patient.ID)
      current_index <- which(
        processed_patient_ids == state$current_processed_patient_id
      )
      next_index <- ifelse(
        current_index == length(processed_patient_ids),
        1, current_index + 1
      )
      state$current_processed_patient_id <- processed_patient_ids[next_index]
      output$processed_id_selector <- renderUI({
        selectizeInput("processed_id_select",
          "Select Patient ID",
          choices = unique(state$preprocessed_data[[1]]$Patient.ID),
          selected = state$current_processed_patient_id
        )
      })
    })

    observeEvent(input$previous_processed_patient, {
      processed_patient_ids <- unique(state$preprocessed_data[[1]]$Patient.ID)
      current_index <- which(
        processed_patient_ids == state$current_processed_patient_id
      )
      prev_index <- ifelse(
        current_index == 1,
        length(processed_patient_ids),
        current_index - 1
      )
      state$current_processed_patient_id <- processed_patient_ids[prev_index]
      output$processed_id_selector <- renderUI({
        selectizeInput("processed_id_select",
          "Select Patient ID",
          choices = unique(state$preprocessed_data[[1]]$Patient.ID),
          selected = state$current_processed_patient_id
        )
      })
    })

    observeEvent(input$preprocess_data, {
      if (is.null(state$data_set)) {
        showNotification(
          "Cannot process data without a dataset",
          type = "message",
          duration = 3
        )
        return()
      }
      showModal(modalDialog(
        title = "Processing Data",
        div(
          "Processing data, this step may take up to a few minutes.",
          div(
            withLoader(
              img(src = "https://i.gifer.com/7efs.gif", height = "200px"),
              type = "html"
            ),
            style = paste(
              "max-height: 150px; overflow: hidden;",
              "display: flex; justify-content: center;",
              "padding-top: 30px;"
            )
          )
        ),
        footer = NULL,
        easyClose = FALSE,
        size = "m"
      ))
      state$preprocessed_data <- prepareDataset(state$data_set, input)
      if (is.null(state$preprocessed_data)) {
        showNotification(
          "Dataprocessing returned no valid rows",
          type = "message",
          duration = 3
        )
        removeModal()
        return()
      }
      if (state$preprocessed_data[[4]] != 0) {
        showNotification(
          paste0(
            "Removed ",
            as.character(state$preprocessed_data[[4]]),
            " line(s) from Dataset due to invalid genomic positions."
          ),
          type = "message",
          duration = 10
        )
      }
      if (state$preprocessed_data[[5]] != 0) {
        showNotification(
          paste0(
            "Removed ",
            as.character(state$preprocessed_data[[5]]),
            " line(s) due to undescribed genomic regions."
          ),
          type = "message",
          duration = 10
        )
      }
      if (state$preprocessed_data[[8]] != 0) {
        showNotification(
          paste0(
            "Removed ",
            as.character(state$preprocessed_data[[8]]),
            " patient(s) from Dataset due to invalid tree structures."
          ),
          type = "message",
          duration = 10
        )
      }
      state$current_processed_patient_id <- unique(
        state$preprocessed_data[[1]]$Patient.ID
      )[[1]]
      output$processed_id_selector <- renderUI({
        selectizeInput("processed_id_select",
          "Select Patient ID",
          choices = unique(state$preprocessed_data[[1]]$Patient.ID),
          selected = state$current_processed_patient_id
        )
      })
      set_data_statistics_panel_to_blank(output, state)
      set_signature_computation_panel_to_blank(output, state)
      set_trajectory_approximation_panel_to_blank(output, state)
      removeModal()
      generate_processed_data_main_panel(state, output)
    })

    ######################################
    # Eventhandlers in the data statistics tab
    ######################################

    observeEvent(input$generate_stats, {
      if (is.null(state$data_set)) {
        showNotification(
          "Cannot generate Statistics without a dataset.",
          type = "message",
          duration = 3
        )
      } else {
        if (is.null(state$preprocessed_data)) {
          showNotification(
            "Cannot generate Statistics without a processed dataset.",
            type = "message",
            duration = 3
          )
          return()
        }
        stats <- generateStatistics(state$preprocessed_data[[2]])
        output$clone_count_plot <- renderPlot({
          clone_data <- as.data.frame(stats$UniqueCloneCountFrequency)
          colnames(clone_data) <- c("CloneCount", "DatasetCount")
          ggplot(
            clone_data,
            aes(x = .data$CloneCount, y = .data$DatasetCount)
          ) +
            geom_col(fill = "#00d8c6") +
            labs(x = "Number of Unique Clones", y = "Number of Patients") +
            theme_bw() +
            theme(
              plot.title = element_text(size = 18, face = "bold"),
              axis.title.x = element_text(size = 16, face = "bold"),
              axis.title.y = element_text(size = 16, face = "bold"),
              axis.text = element_text(size = 12),
              legend.title = element_text(size = 16, face = "bold"),
              legend.text = element_text(size = 16)
            )
        })

        output$variant_count_plot <- renderPlot({
          variant_data <- as.data.frame(stats$RowCountFrequency)
          colnames(variant_data) <- c("VariantCount", "DatasetCount")
          ggplot(
            variant_data,
            aes(x = .data$VariantCount, y = .data$DatasetCount)
          ) +
            geom_col(fill = "darkorange") +
            labs(x = "Number of Variants", y = "Number of Patients") +
            theme_bw() +
            theme(
              plot.title = element_text(size = 18, face = "bold"),
              axis.title.x = element_text(size = 16, face = "bold"),
              axis.title.y = element_text(size = 16, face = "bold"),
              axis.text = element_text(size = 12),
              legend.title = element_text(size = 16, face = "bold"),
              legend.text = element_text(size = 16)
            )
        })

        output$variants_multi_occurrence_plot <- renderPlot({
          variant_freq <- table(state$preprocessed_data[[2]]$Variant)
          freq_of_freqs <- table(variant_freq)
          par(mar = c(5, 4, 1, 1))
          layout(matrix(c(1), nrow = 1, ncol = 1))
          plot.new()
          barplot(freq_of_freqs,
            col = "darkorange",
            border = "black",
            las = 1,
            cex.axis = 1.2,
            font.axis = 2
          )
          title(
            xlab = "Number of Occurrences",
            ylab = "Number of Variants",
            font.lab = 2
          )
        })

        chromosomeData <- stats$ChromosomePercentage
        if (!is.null(chromosomeData)) {
          output$chromosome_distribution_plot <- renderPlot({
            chromosome_order <- as.character(c(seq_len(22), "X", "Y"))
            sortedData <- chromosomeData[chromosome_order]
            par(mar = c(5, 4, 1, 1))
            layout(matrix(c(1), nrow = 1, ncol = 1))
            barplot(sortedData,
              names.arg = names(sortedData),
              col = "#00d8c6",
              border = "black",
              las = 1,
              cex.axis = 1.2,
              font.axis = 2
            )
            title(
              xlab = "Chromosome",
              ylab = "Percentage of Variants",
              font.lab = 2
            )
          })
        }

        output$CNV_SNV_ratio_plot <- renderUI({
          snv_perc <- stats$SNVINDELPercentage
          cnv_perc <- stats$CNVPercentage
          na_perc <- stats$NAPercentage

          snv_width <- paste0(snv_perc, "%")
          cnv_width <- paste0(cnv_perc, "%")
          na_width <- paste0(na_perc, "%")

          div(
            style = "width: 100%;",
            div(
              style = "height: 30px; display: flex; border: 1px solid black;",
              div(
                style = paste0(
                  "width:", snv_width,
                  "; background-color: darkorange;"
                )
              ),
              div(
                style = paste0(
                  "width:", na_width,
                  "; background-color: grey;"
                )
              ),
              div(
                style = paste0(
                  "width:", cnv_width,
                  "; background-color: #00d8c6;"
                )
              )
            ),
            div(
              style = "display: flex; font-weight: bold; margin-top: 5px;",
              div(
                style = paste(
                  "width: 33.33%; text-align: left;",
                  "padding-left: 10px; color: darkorange;"
                ),
                paste("SNVs/INDELs - ", round(snv_perc, 2), "%")
              ),
              div(
                style = "width: 33.33%; text-align: center; color: grey;",
                paste("NA - ", round(na_perc, 2), "%")
              ),
              div(
                style = paste(
                  "width: 33.33%; text-align: right;",
                  "padding-right: 10px; color: #00d8c6;"
                ),
                paste("CNVs - ", round(cnv_perc, 2), "%")
              )
            )
          )
        })

        output$top_variants_plot <- renderPlot({
          topVariants <- stats$Top10Variants
          topVariantsPercentage <- stats$Top10VariantsPercentage
          topVariantsData <- data.frame(
            Variant = names(topVariants),
            Frequency = as.numeric(topVariants),
            Percentage = as.numeric(topVariantsPercentage)
          )
          ggplot(
            topVariantsData,
            aes(
              x = .data$Frequency,
              y = reorder(.data$Variant, .data$Frequency)
            )
          ) +
            geom_bar(stat = "identity", fill = "#00d8c6") +
            geom_text(
              aes(label = paste0(round(.data$Percentage, 1), "%")),
              hjust = -0.2,
              size = 4
            ) +
            labs(x = "Frequency", y = "Variants") +
            theme_bw() +
            theme(
              axis.text.y = element_text(size = 12), ,
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.x = element_blank(),
              axis.title = element_text(size = 14, face = "bold"),
              plot.title = element_text(size = 16, face = "bold")
            )
        })
      }
    })

    ######################################
    # Eventhandlers in the signature computation tab
    ######################################

    observeEvent(input$compute_signatures, {
      if (is.null(state$data_set)) {
        showNotification("No dataset loaded yet.
                        Please check the Data Editor tab.",
          type = "message",
          duration = 3
        )
        return()
      }
      if (is.null(state$preprocessed_data)) {
        showNotification("Dataset has not been processed yet.
                        Please check the Processed Data tab.",
          type = "message",
          duration = 3
        )
        return()
      }

      showModal(
        modalDialog(
          title = "Computing Signatures",
          div(
            "Computing Signatures, this process may take a while...",
            div(
              withLoader(
                img(
                  src = "https://i.gifer.com/7efs.gif",
                  height = "200px"
                ),
                type = "html"
              ),
              style = "max-height: 150px;
                      overflow: hidden;
                      display: flex;
                      justify-content: center;
                      padding-top: 30px;"
            )
          ),
          footer = NULL,
          easyClose = FALSE,
          size = "m"
        )
      )

      signature_information <- compute_signatures(
        state$preprocessed_data,
        input
      )
      state$cnv_region_result <- NULL
      output$cnv_region_panel <- NULL

      # Preserve an empty seed state for independent-root reconstruction.
      if (
        is.null(signature_information) ||
          length(signature_information) < 1 ||
          is.null(signature_information[[1]]) ||
          nrow(signature_information[[1]]) == 0
      ) {
        removeModal()
        # Independent root branches can exist without any linear signatures.
        # Keep a typed empty seed table so trajectory computation remains usable.
        state$signature_data_set <- data.frame(
          `Variant 1` = character(), `Number of occurrences` = numeric(),
          Patients = character(), `p value` = numeric(),
          `Relative occurrence among patients with identical variants` = numeric(),
          check.names = FALSE)
        state$filtered_signature_data_set <- NULL
        state$signature_side_information <- if (length(signature_information) >= 2L) signature_information[[2]] else NULL
        state$signature_type <- if (input$radio_multiple_variants_per_clone) 2L else 1L
        showNotification(
          paste("No linear signatures meet the current criteria.",
            "Independent root branches can still be computed in the trajectory tab."),
          type = "message",
          duration = 5
        )
        return()
      }

      if (input$radio_multiple_variants_per_clone) {
        state$signature_type <- 2
      } else {
        state$signature_type <- 1
      }

      state$signature_data_set <- signature_information[[1]]
      state$filtered_signature_data_set <- NULL
      state$signature_side_information <- signature_information[[2]]

      output$signature_computation_main_panel <- renderUI({
        DT::dataTableOutput("signature_table")
      })

      output$signature_table <- DT::renderDataTable({
        datatable(
          state$signature_data_set,
          options = list(
            scrollX = TRUE,
            dom = "tp",
            paging = TRUE,
            pageLength = 20,
            lengthMenu = c(10, 25, 50, 100),
            columnDefs = list(
              list(className = "dt-left", targets = "_all"),
              list(visible = FALSE, targets = "Patients")
            ),
            fixedHeader = TRUE
          ),
          editable = FALSE,
          rownames = FALSE,
          selection = "single"
        )
      })

      set_trajectory_approximation_panel_to_blank(output, state)
      removeModal()
    })

    observeEvent(input$apply_signature_filter, {
      if (is.null(state$filtered_signature_data_set)) {
        state$filtered_signature_data_set <- state$signature_data_set
        state$filtered_signature_data_set <-
          state$filtered_signature_data_set %>%
          mutate(
            occurrences = as.numeric(.data[["Number of occurrences"]])
          ) %>%
          filter(
            !is.na(.data$occurrences) &
              .data$occurrences >= input$absolute_frequency_signatures
          )
        state$filtered_signature_data_set <-
          state$filtered_signature_data_set %>%
          select(-all_of("occurrences"))
        if (input$signature_search != "") {
          variant_columns <- grep("^Variant ", names(state$signature_data_set))
          searchable_signature_columns <- state$filtered_signature_data_set[
            ,
            variant_columns,
            drop = FALSE
          ]
          matching_signatures <- apply(
            searchable_signature_columns,
            1,
            function(row) {
              any(grepl(input$signature_search, row, ignore.case = TRUE))
            }
          )
          state$filtered_signature_data_set <-
            state$filtered_signature_data_set[matching_signatures, ]
        }
        if ("p value" %in% names(state$filtered_signature_data_set)) {
          state$filtered_signature_data_set <-
            state$filtered_signature_data_set %>%
            filter(.data[["p value"]] < input$min_p_value)
          state$filtered_signature_data_set <-
            state$filtered_signature_data_set %>%
            filter(
              .data[[
                "Relative occurrence among patients with identical variants"
              ]] >=
                input$relative_frequency
            )
        }
        output$signature_computation_main_panel <- renderUI({
          DT::dataTableOutput("signature_table")
        })
        output$signature_table <- DT::renderDataTable({
          dt <- datatable(state$filtered_signature_data_set,
            options = list(
              scrollX = TRUE,
              dom = "tp",
              paging = TRUE,
              pageLength = 20,
              lengthMenu = c(10, 25, 50, 100),
              columnDefs = list(
                list(className = "dt-left", targets = "_all"),
                list(visible = FALSE, targets = "Patients")
              ),
              fixedHeader = TRUE
            ),
            editable = FALSE,
            rownames = FALSE,
            selection = "single"
          )
        })
        updateActionButton(
          session,
          "apply_signature_filter",
          label = "Disable Filter"
        )
      } else {
        updateActionButton(
          session,
          "apply_signature_filter",
          label = "Enable Filter"
        )
        state$filtered_signature_data_set <- NULL
        output$signature_computation_main_panel <- renderUI({
          DT::dataTableOutput("signature_table")
        })
        output$signature_table <- DT::renderDataTable({
          dt <- datatable(state$signature_data_set,
            options = list(
              scrollX = TRUE,
              dom = "tp",
              paging = TRUE,
              pageLength = 20,
              lengthMenu = c(10, 25, 50, 100),
              columnDefs = list(
                list(className = "dt-left", targets = "_all"),
                list(visible = FALSE, targets = "Patients")
              ),
              fixedHeader = TRUE
            ),
            editable = FALSE,
            rownames = FALSE,
            selection = "single"
          )
        })
      }
    })

    observeEvent(input$signature_length_setter, {
      if (is.numeric(input$signature_length_setter)) {
        if (input$signature_length_setter < 2) {
          updateNumericInput(
            session,
            inputId = "signature_length_setter",
            value = 2
          )
        }
      } else {
        updateNumericInput(
          session,
          inputId = "signature_length_setter",
          value = 2
        )
      }
    })

    observeEvent(input$postprocess_cnv_regions, {
      if (is.null(state$signature_side_information)) {
        showNotification("Compute signatures before analyzing CNV regions.", type = "message")
        return()
      }
      result <- tryCatch(postprocess_cnv_regions(
        list(state$signature_data_set, state$signature_side_information),
        minimum_support = input$minimum_support_of_signatures),
        error = function(error) { showNotification(conditionMessage(error), type = "error", duration = 10); NULL })
      if (is.null(result)) return()
      state$cnv_region_result <- result
      output$cnv_region_panel <- renderUI(tagList(
        h3("Shared CNV region signatures"),
        p("Each region is an alternative description of one CNV event. Counts are unique patients; significance must be assessed separately."),
        DT::dataTableOutput("cnv_region_table")))
      output$cnv_region_table <- DT::renderDataTable(DT::datatable(
        state$cnv_region_result$signatures, rownames = FALSE,
        options = list(scrollX = TRUE, pageLength = 20)))
    })
    output$export_cnv_regions <- downloadHandler(
      filename = function() "cnv-region-signatures.csv",
      content = function(file) {
        req(state$cnv_region_result)
        write.csv(state$cnv_region_result$signatures, file, row.names = FALSE)
      })

    observeEvent(input$estimate_significance, {
      if (
        is.null(state$signature_data_set) &&
          is.null(state$filtered_signature_data_set)
      ) {
        showNotification("No signatures computed yet.",
          type = "message",
          duration = 3
        )
        return()
      }
      showModal(modalDialog(
        title = "Computing Significance",
        div(
          "Estimating significance, this process may take a while...",
          div(
            withLoader(
              img(src = "https://i.gifer.com/7efs.gif", height = "200px"),
              type = "html"
            ),
            style = paste(
              "max-height: 150px; overflow: hidden;",
              "display: flex; justify-content: center;",
              "padding-top: 30px;"
            )
          )
        ),
        footer = NULL,
        easyClose = FALSE,
        size = "m"
      ))
      stats <- compute_alt_statistics(
        state$signature_data_set,
        state$signature_side_information[[1]],
        state$signature_side_information[[2]],
        state$signature_side_information[[3]],
        input
      )
      state$signature_data_set <- add_stats_to_dataframe(
        state$signature_data_set,
        stats
      )
      state$signature_data_set <- state$signature_data_set %>%
        mutate(`p value` = as.numeric(.data[["p value"]]))
      removeModal()
    })

    observeEvent(input$signature_table_rows_selected, {
      selected_row <- input$signature_table_rows_selected
      true_row_index <- input$signature_table_rows_all[selected_row]
      if (length(selected_row) == 0) {
        return()
      }
      row_data <- state$signature_data_set[true_row_index, , drop = FALSE]
      variant_cols <- grep("^Variant \\d+$", names(row_data), value = TRUE)
      variant_values <- unlist(row_data[variant_cols])
      processed_df <- state$preprocessed_data[[1]]

      patients_string <- row_data$Patients
      patients_string <- gsub('c\\(|\\)|"', "", row_data$Patients)
      affected_patients <- strsplit(patients_string, ",")[[1]]
      affected_patients <- trimws(affected_patients)

      unique_map <- !duplicated(variant_values)
      variant_values <- variant_values[unique_map]

      bar_data <- data.frame(
        Variant = variant_values,
        Value = as.character(variant_values),
        Percent = vapply(variant_values, function(val) {
          percent <- mean(processed_df$Variant == val, na.rm = TRUE) * 100
          return(percent)
        }, numeric(1))
      )

      bar_data <- na.omit(bar_data)

      output$variant_barplot <- renderPlot({
        ggplot(
          bar_data,
          aes(x = .data$Variant, y = .data$Percent, fill = .data$Variant)
        ) +
          geom_bar(stat = "identity") +
          scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
          theme_minimal(base_size = 18) +
          theme(
            axis.title = element_text(face = "bold"),
            axis.text = element_text(face = "bold", size = 16),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.title = element_text(face = "bold"),
            legend.text = element_text(size = 14),
            legend.position = "none"
          )
      })

      affected_patients_str <- paste(affected_patients, collapse = ", ")

      showModal(modalDialog(
        title = tags$div(style = "font-weight: bold;", "Signature Statistics"),
        div(
          h3(
            "Affected Patients",
            style = "text-align: left; font-weight: bold; font-size: 22px;"
          ),
          p(affected_patients_str),
          hr(
            style = paste(
              "border: 1px solid #D3D3D3; margin-top: 20px;",
              "margin-bottom: 20px;"
            )
          ),
          h3(
            "Variant Frequencies",
            style = paste(
              "text-align: left; font-weight: bold;",
              "font-size: 22px; margin-top: 20px;"
            )
          )
        ),
        plotOutput("variant_barplot"),
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
    })

    ######################################
    # Eventhandlers in the trajectory computation tab
    ######################################

    observeEvent(input$compute_trajectories, {
      if (
        !is.null(state$signature_data_set) &&
          "p value" %in% colnames(state$signature_data_set)
      ) {
        showModal(modalDialog(
          title = "Computing Trajectories",
          div(
            paste(
              "The signatures defined through your filters will now be",
              "recombined into full branching trajectories and analyzed.",
              "For progress indication check the console"
            ),
            div(
              withLoader(
                img(src = "https://i.gifer.com/7efs.gif", height = "200px"),
                type = "html"
              ),
              style = paste(
                "max-height: 150px; overflow: hidden;",
                "display: flex; justify-content: center;",
                "padding-top: 30px;"
              )
            )
          ),
          footer = NULL,
          easyClose = FALSE,
          size = "m"
        ))

        filtered_data <- state$signature_data_set
        filtered_data <- filtered_data %>%
          mutate(
            occurrences = as.numeric(.data[["Number of occurrences"]])
          ) %>%
          filter(
            !is.na(.data$occurrences) &
              .data$occurrences >= input$absolute_frequency_signatures
          )
        filtered_data <- filtered_data %>% select(-all_of("occurrences"))

        if (input$signature_search != "") {
          variant_columns <- grep("^Variant ", names(state$signature_data_set))
          filtered_data <- filtered_data[
            apply(
              filtered_data[, variant_columns, drop = FALSE],
              1,
              function(row) {
                any(grepl(input$signature_search, row, ignore.case = TRUE))
              }
            ),
          ]
        }

        if ("p value" %in% names(filtered_data)) {
          filtered_data <- filtered_data %>%
            filter(.data[["p value"]] < input$min_p_value)
          filtered_data <- filtered_data %>%
            filter(
              .data[[
                "Relative occurrence among patients with identical variants"
              ]] >=
                input$relative_frequency
            )
        }
        trajectory_result <- compute_trajectories(
          filtered_data,
          state$preprocessed_data,
          input,
          state$signature_type
        )

        if (is.null(trajectory_result) || nrow(trajectory_result) == 0) {
          removeModal()
          showNotification(
            "No trajectories meet the current criteria.",
            type = "message",
            duration = 5
          )
          return()
        }

        state$trajectory_data_set <- trajectory_result
        state$filtered_trajectory_data_set <- NULL

        removeModal()

        output$trajectory_approximation_main_panel <- renderUI({
          DT::dataTableOutput("trajectory_table")
        })
        output$trajectory_table <- DT::renderDataTable({
          datatable(
            state$trajectory_data_set,
            options = list(
              scrollX = TRUE,
              dom = "tp",
              paging = TRUE,
              pageLength = 20,
              lengthMenu = c(10, 25, 50, 100),
              columnDefs = list(
                list(className = "dt-left", targets = "_all"),
                list(visible = FALSE, targets = c("Patients"))
              ),
              fixedHeader = TRUE
            ),
            editable = FALSE,
            rownames = FALSE,
            selection = "none"
          )
        })
      } else {
        showNotification(
          paste(
            "Cannot process trajectories without precomputed and",
            "estimated signatures"
          ),
          type = "message",
          duration = 3
        )
      }
    })

    observeEvent(input$recombination_number_setter, {
      if (is.numeric(input$recombination_number_setter)) {
        if (input$recombination_number_setter < 2) {
          updateNumericInput(
            session,
            inputId = "recombination_number_setter",
            value = 1
          )
        }
      } else {
        updateNumericInput(
          session,
          inputId = "recombination_number_setter",
          value = 1
        )
      }
    })

    observeEvent(input$estimate_significant_trajectories, {
      if (is.null(state$trajectory_data_set)) {
        showNotification("No trajectories computed yet.",
          type = "message",
          duration = 3
        )
        return()
      } else {
        showModal(modalDialog(
          title = "Computing trajectory significance measures",
          div(
            paste(
              "The trajectories computed will now be assessed through",
              "the settings configured."
            ),
            div(
              withLoader(
                img(src = "https://i.gifer.com/7efs.gif", height = "200px"),
                type = "html"
              ),
              style = paste(
                "max-height: 150px; overflow: hidden;",
                "display: flex; justify-content: center;",
                "padding-top: 30px;"
              )
            )
          ),
          footer = NULL,
          easyClose = FALSE,
          size = "m"
        ))
        state$trajectory_data_set <- estimate_significance(
          state$trajectory_data_set,
          state$preprocessed_data[[1]],
          input,
          state$signature_type
        )
        removeModal()
      }
    })

    observeEvent(input$apply_trajectory_filter, {
      if (is.null(state$filtered_trajectory_data_set)) {
        state$filtered_trajectory_data_set <- state$trajectory_data_set
        state$filtered_trajectory_data_set <-
          state$filtered_trajectory_data_set %>%
          mutate(
            occurrences = as.numeric(
              str_extract(
                str_trim(as.character(.data[["Number of occurrences"]])),
                "^\\d+"
              )
            )
          ) %>%
          filter(
            !is.na(.data$occurrences) &
              .data$occurrences >= input$absolute_frequency_trajectories
          )
        state$filtered_trajectory_data_set <-
          state$filtered_trajectory_data_set %>%
          select(-all_of("occurrences"))
        if (input$variant_search_trajectories != "") {
          state$filtered_trajectory_data_set <-
            state$filtered_trajectory_data_set %>%
            dplyr::group_by(.data$ID) %>%
            dplyr::filter(
              any(grepl(
                input$variant_search_trajectories,
                .data$Variant,
                ignore.case = TRUE
              ))
            ) %>%
            dplyr::ungroup()
        }
        if ("p value" %in% names(state$filtered_trajectory_data_set)) {
          state$filtered_trajectory_data_set <-
            state$filtered_trajectory_data_set %>%
            filter(.data[["p value"]] < input$min_p_value_trajectories)
          state$filtered_trajectory_data_set <-
            state$filtered_trajectory_data_set %>%
            filter(
              .data[[
                "Relative occurrence among patients with identical variants"
              ]] >=
                input$relative_frequency_trajectories
            )
        }
        output$trajectory_approximation_main_panel <- renderUI({
          DT::dataTableOutput("trajectory_table")
        })
        output$trajectory_table <- DT::renderDataTable({
          dt <- datatable(state$filtered_trajectory_data_set,
            options = list(
              scrollX = TRUE,
              dom = "tp",
              paging = TRUE,
              pageLength = 20,
              lengthMenu = c(10, 25, 50, 100),
              columnDefs = list(
                list(className = "dt-left", targets = "_all"),
                list(visible = FALSE, targets = c("Patients", "Shape"))
              ),
              fixedHeader = TRUE
            ),
            editable = FALSE,
            rownames = FALSE,
            selection = "single"
          )
        })
        updateActionButton(
          session,
          "apply_trajectory_filter",
          label = "Disable Filter"
        )
      } else {
        updateActionButton(
          session,
          "apply_trajectory_filter",
          label = "Enable Filter"
        )
        output$trajectory_approximation_main_panel <- renderUI({
          DT::dataTableOutput("trajectory_table")
        })
        output$trajectory_table <- DT::renderDataTable({
          datatable(
            state$trajectory_data_set,
            options = list(
              scrollX = TRUE,
              dom = "tp",
              paging = TRUE,
              pageLength = 20,
              lengthMenu = c(10, 25, 50, 100),
              columnDefs = list(
                list(className = "dt-left", targets = "_all"),
                list(visible = FALSE, targets = c("Patients", "Shape"))
              ),
              fixedHeader = TRUE
            ),
            editable = FALSE,
            rownames = FALSE,
            selection = "none"
          )
        })
        state$filtered_trajectory_data_set <- NULL
      }
    })

    observeEvent(input$export_trajectories, {
      if (is.null(state$filtered_trajectory_data_set)) {
        write.csv(
          state$trajectory_data_set,
          file = "clevRact_trajectories.csv",
          row.names = FALSE
        )
      } else {
        write.csv(
          state$filtered_trajectory_data_set,
          file = "clevRact_trajectories.csv",
          row.names = FALSE
        )
      }
      message(
        "Data frame successfully exported to ",
        "clevRact_trajectories.csv"
      )
    })

    ######################################
    # Universal handlers
    ######################################

    observeEvent(input$file, {
      req(input$file)
      removeModal()
      file_path <- input$file$datapath
      if (handle_loaded_dataframe(state, output, input)) {
        set_processed_data_panel_to_blank(output, state)
        set_data_statistics_panel_to_blank(output, state)
        set_signature_computation_panel_to_blank(output, state)
        set_trajectory_approximation_panel_to_blank(output, state)
      }
    })

    output$downloadData <- downloadHandler(
      filename = function() {
        paste("my_data_", Sys.Date(), Sys.time(), ".csv")
      },
      content = function(file) {
        write.table(
          state$data_set,
          file = file,
          row.names = FALSE,
          dec = ".",
          sep = ";",
          quote = FALSE
        )
      }
    )

    ######################################
    # Eventhandlers for state variables
    ######################################

    observeEvent(state$current_patient_id, {
      if (!is.null(state$current_patient_id)) {
        generate_editor_main_panel(state, input, output)
      }
    })

    observeEvent(state$current_processed_patient_id, {
      if (!is.null(state$current_processed_patient_id)) {
        generate_processed_data_main_panel(state, output)
      }
    })
  }
}
