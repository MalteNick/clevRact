tutorial_input <- function(exact = FALSE, within = FALSE, roots = TRUE) {
  list(dropdown_variant_generation = "Yes, using GRCh37",
    dropdown_summary_snvs = if (exact) "Use exact SNVs" else "Summarize on gene level",
    dropdown_summary_cnvs = "Use exact CNVs", dropdown_cnv_clone_summary = FALSE,
    dropdown_handler_reoccurrence = "No differentiation",
    minimum_support_of_signatures = 2L, signature_length_setter = 2L,
    radio_multiple_variants_per_clone = within,
    compute_independent_root_branches = roots,
    multiple_identical_variants_for_branches = TRUE,
    minimum_support_of_trajectories = 2L, recombination_number_setter = 2L,
    significance_estimation_method = "Sandmann",
    p_value_adjustment_trajectories = "None")
}
tutorial_data <- function() {
  read.csv(system.file("extdata", "synthetic_10_patients.csv", package = "clevRact"),
    sep = ";")
}
tutorial_signature_support <- function(signatures, first, second) {
  rows <- gsub("\\s+", "", signatures$`Variant 1`) == first &
    gsub("\\s+", "", signatures$`Variant 2`) == second
  as.numeric(signatures$`Number of occurrences`[rows])
}
tutorial_trajectory_support <- function(result, labels, parents) {
  target <- data.frame(Variant = labels, Clone = seq_along(labels), Parent = parents)
  result$Variant <- gsub("\\s+", "", result$Variant)
  groups <- split(result, result$ID)
  matches <- Filter(function(x) {
    nrow(x) == nrow(target) && compare_trajectories(x, target)
  }, groups)
  vapply(matches, function(x) as.numeric(x$`Number of occurrences`[1]), numeric(1))
}

test_that("the tutorial includes ten valid histories and both time points", {
  data <- tutorial_data()
  expect_equal(nrow(data), 46L)
  expect_equal(length(unique(data$Patient.ID)), 10L)
  expect_true(all(c("T1", "T2") %in% names(data)))
  expect_true(all(data$T2))
  expect_true(any(!data$T1))
  prepared <- prepareDataset(data, tutorial_input())
  expect_equal(nrow(prepared[[1]]), nrow(data))
  expect_equal(prepared$noi_trees, 0L)
  for (gene in c("PHF6", "DDX3Y")) {
    rows <- prepared[[1]][prepared[[1]]$Variant == gene, ]
    expect_equal(length(unique(rows$Patient.ID)), 2L)
  }
})

test_that("gene and exact-position encoding give the documented linear support", {
  for (exact in c(FALSE, TRUE)) {
    input <- tutorial_input(exact = exact)
    prepared <- prepareDataset(tutorial_data(), input)
    signatures <- compute_signatures(prepared, input)
    a <- if (exact) "9_139390656" else "NOTCH1"
    b <- if (exact) "9_139397640" else "NOTCH1"
    f <- if (exact) "4_153245299" else "FBXW7"
    expect_equal(tutorial_signature_support(signatures[[1]], a, f), 10)
    expect_equal(tutorial_signature_support(signatures[[1]], a, b), 2)
    input$radio_multiple_variants_per_clone <- TRUE
    relaxed <- compute_signatures(prepared, input)
    expect_equal(tutorial_signature_support(relaxed[[1]], a, b), 4)
    values <- list()
    for (method in c("Sandmann", "Sandmann Nick")) {
      input$significance_estimation_method <- method
      stats <- compute_alt_statistics(signatures[[1]], prepared[[1]],
        signatures[[2]][[2]], signatures[[2]][[3]], input)
      scored <- add_stats_to_dataframe(signatures[[1]], stats)
      values[[method]] <- as.numeric(scored$`p value`)
      expect_true(all(is.finite(values[[method]])))
      expect_true(all(values[[method]] >= 0 & values[[method]] <= 1))
    }
    expect_equal(isTRUE(all.equal(values[[1]], values[[2]])), exact)
  }
})

test_that("tutorial CNV regions pool singleton intervals and retain event types", {
  input <- tutorial_input()
  prepared <- prepareDataset(tutorial_data(), input)
  signatures <- compute_signatures(prepared, input)
  expect_equal(tutorial_signature_support(signatures[[1]],
    "CNV_del_9_21600000_22050000", "NOTCH1"), 2)
  result <- postprocess_cnv_regions(signatures, minimum_support = 2L)
  expect_equal(tutorial_signature_support(result$signatures,
    "CNV_del_region_9_21800000_22050000", "NOTCH1"), 8)
  for (type in c("CNV_dup", "CNV_loh")) {
    expect_equal(tutorial_signature_support(result$signatures,
      paste0(type, "_region_9_24000000_25000000"), "NOTCH1"), 2)
  }
  expect_false("p value" %in% names(result$signatures))
})

test_that("tutorial branching distinguishes independent and connected structures", {
  input <- tutorial_input()
  prepared <- prepareDataset(tutorial_data(), input)
  signatures <- compute_signatures(prepared, input)
  stats <- compute_alt_statistics(signatures[[1]], prepared[[1]],
    signatures[[2]][[2]], signatures[[2]][[3]], input)
  signatures[[1]] <- add_stats_to_dataframe(signatures[[1]], stats)
  result <- compute_trajectories(signatures[[1]], prepared, input, 1L)
  expect_equal(unname(tutorial_trajectory_support(result,
    c("NOTCH1", "NOTCH1"), c(0, 0))), 6)
  for (adjustment in c("None", "Benjamini-Hochberg")) {
    input$p_value_adjustment_trajectories <- adjustment
    scored <- estimate_significance(result, prepared[[1]], input, 1L)
    expect_true(all(is.finite(scored$`p value`)))
    expect_true(all(scored$`p value` >= 0 & scored$`p value` <= 1))
  }
  input$compute_independent_root_branches <- FALSE
  connected <- compute_trajectories(signatures[[1]], prepared, input, 1L)
  expect_length(tutorial_trajectory_support(connected,
    c("NOTCH1", "NOTCH1"), c(0, 0)), 0L)
  expect_equal(unname(tutorial_trajectory_support(connected,
    c("CNV_del_9_21600000_22050000", "NOTCH1", "NOTCH1", "FBXW7"),
    c(0, 1, 1, 2))), 2)
  input$compute_independent_root_branches <- TRUE
  input$dropdown_summary_snvs <- "Use exact SNVs"
  prepared <- prepareDataset(tutorial_data(), input)
  signatures <- compute_signatures(prepared, input)
  stats <- compute_alt_statistics(signatures[[1]], prepared[[1]],
    signatures[[2]][[2]], signatures[[2]][[3]], input)
  signatures[[1]] <- add_stats_to_dataframe(signatures[[1]], stats)
  result <- compute_trajectories(signatures[[1]], prepared, input, 1L)
  expect_equal(unname(tutorial_trajectory_support(result,
    c("9_139390656", "9_139397640"), c(0, 0))), 4)
  expect_equal(unname(tutorial_trajectory_support(result,
    c("9_139390656", "9_139399365"), c(0, 0))), 2)
})

test_that("the synthetic cohort can be analysed through the GUI controls", {
  shiny::testServer(.build_server(), {
    path <- system.file("extdata", "synthetic_10_patients.csv", package = "clevRact")
    session$setInputs(toggle_view = FALSE,
      file = data.frame(name = basename(path), type = "text/csv",
        size = file.info(path)$size, datapath = path))
    expect_equal(length(unique(state$data_set$Patient.ID)), 10L)
    session$setInputs(dropdown_variant_generation = "Yes, using GRCh37",
      dropdown_summary_snvs = "Summarize on gene level",
      dropdown_summary_cnvs = "Use exact CNVs", dropdown_cnv_clone_summary = "FALSE",
      dropdown_handler_reoccurrence = "No differentiation", preprocess_data = 1L)
    expect_equal(nrow(state$preprocessed_data[[1]]), 46L)
    session$setInputs(minimum_support_of_signatures = 2L,
      signature_length_setter = 2L, radio_multiple_variants_per_clone = FALSE,
      compute_signatures = 1L)
    expect_equal(tutorial_signature_support(state$signature_data_set,
      "NOTCH1", "FBXW7"), 10)
    for (method in c("Sandmann", "Sandmann Nick")) {
      session$setInputs(significance_estimation_method = method,
        estimate_significance = if (method == "Sandmann") 1L else 2L)
      expect_true(all(is.finite(state$signature_data_set$`p value`)))
    }
    session$setInputs(significance_estimation_method = "Sandmann",
      estimate_significance = 3L, postprocess_cnv_regions = 1L)
    expect_equal(tutorial_signature_support(state$cnv_region_result$signatures,
      "CNV_del_region_9_21800000_22050000", "NOTCH1"), 8)
    session$setInputs(absolute_frequency_signatures = 2,
      signature_search = "", min_p_value = 1, relative_frequency = 0,
      minimum_support_of_trajectories = 2L, recombination_number_setter = 2L,
      compute_independent_root_branches = "TRUE",
      multiple_identical_variants_for_branches = "TRUE", compute_trajectories = 1L)
    expect_equal(unname(tutorial_trajectory_support(state$trajectory_data_set,
      c("NOTCH1", "NOTCH1"), c(0, 0))), 6)
    session$setInputs(p_value_adjustment_trajectories = "Benjamini-Hochberg",
      estimate_significant_trajectories = 1L)
    expect_true(all(is.finite(state$trajectory_data_set$`p value`)))
    session$setInputs(compute_independent_root_branches = "FALSE",
      compute_trajectories = 2L)
    expect_length(tutorial_trajectory_support(state$trajectory_data_set,
      c("NOTCH1", "NOTCH1"), c(0, 0)), 0L)
    session$setInputs(radio_multiple_variants_per_clone = TRUE,
      compute_signatures = 2L)
    expect_equal(tutorial_signature_support(state$signature_data_set,
      "NOTCH1", "NOTCH1"), 4)
    session$setInputs(dropdown_summary_snvs = "Use exact SNVs", preprocess_data = 2L)
    session$setInputs(radio_multiple_variants_per_clone = FALSE, compute_signatures = 3L)
    expect_equal(tutorial_signature_support(state$signature_data_set,
      "9_139390656", "4_153245299"), 10)
  })
})
