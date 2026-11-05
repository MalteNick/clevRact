test_that("example input CSV files are installed", {
  examples <- c(
    "synthetic_10_patients.csv"
  )

  paths <- system.file("extdata", examples, package = "clevRact")

  expect_true(all(nzchar(paths)))
  expect_true(all(file.exists(paths)))
})

test_that("example input CSV files contain the required columns", {
  required_columns <- c(
    "Patient.ID",
    "Variant",
    "Clone",
    "Parent",
    "Variant.Type",
    "Chromosome",
    "Start",
    "End"
  )

  examples <- c(
    "synthetic_10_patients.csv"
  )

  for (example in examples) {
    path <- system.file("extdata", example, package = "clevRact")
    data <- utils::read.csv(path, sep = ";", nrows = 1)

    expect_true(
      all(required_columns %in% names(data)),
      info = example
    )
  }
})

test_that("annotation resources are installed", {
  annotation_files <- c(
    "GRCh37.csv",
    "GRCh38.csv",
    "Cytobands_GRCh37.csv",
    "Cytobands_GRCh38.csv",
    "arm_positions.csv"
  )

  paths <- system.file(
    "extdata",
    "annotations",
    annotation_files,
    package = "clevRact"
  )

  expect_true(all(nzchar(paths)))
  expect_true(all(file.exists(paths)))
})
