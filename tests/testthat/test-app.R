test_that("clevRactApp returns a Shiny app object", {
  app <- clevRactApp()

  expect_s3_class(app, "shiny.appobj")
  expect_true(is.function(app$serverFunc))
})

test_that("the GUI exposes optional analysis controls without an About tab", {
  ui <- as.character(.build_ui())
  expect_false(grepl('data-value="About"', ui, fixed = TRUE))
  for (tab in c("Input Data", "Processed Data", "Data Statistics",
    "Linear Trajectories", "Trajectory Reconstruction")) {
    expect_true(grepl(paste0('data-value="', tab, '"'), ui, fixed = TRUE))
  }
  for (id in c("compute_independent_root_branches", "postprocess_cnv_regions",
    "export_cnv_regions", "cnv_region_panel")) {
    expect_true(grepl(paste0('id="', id, '"'), ui, fixed = TRUE))
  }
})
