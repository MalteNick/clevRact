#' Build the clevRact Shiny application
#'
#' Creates the clevRact user interface and server without launching the
#' application. Pass the returned object to [shiny::runApp()] to start it.
#'
#' @return A Shiny application object.
#' @export
#'
#' @examples
#' app <- clevRactApp()
#' if (interactive()) {
#'   shiny::runApp(app)
#' }
clevRactApp <- function() {
  assets <- system.file("www", package = "clevRact")
  if (!nzchar(assets)) {
    stop("The clevRact web assets could not be found.", call. = FALSE)
  }

  shiny::addResourcePath("clevRact-assets", assets)

  shiny::shinyApp(
    ui = .build_ui(),
    server = .build_server()
  )
}
