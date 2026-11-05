.annotation_file <- function(filename) {
  path <- system.file(
    "extdata",
    "annotations",
    filename,
    package = "clevRact"
  )

  if (!nzchar(path)) {
    stop(
      sprintf(
        "The packaged annotation file '%s' could not be found.",
        filename
      ),
      call. = FALSE
    )
  }

  path
}
