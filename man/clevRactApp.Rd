\name{clevRactApp}
\alias{clevRactApp}
\title{Build the clevRact Shiny Application}
\usage{
clevRactApp()
}
\value{
A \code{shiny.appobj} containing the clevRact user interface and server
function.
}
\description{
Creates the clevRact Shiny application for interactive preparation,
visualization, and analysis of longitudinal cancer clonal evolution data.
}
\details{
\code{clevRactApp()} constructs and returns the application without launching
it. This allows callers to choose how the application is run and to configure
options supported by \code{\link[shiny:runApp]{shiny::runApp()}}.

The application accepts semicolon-separated CSV input describing patients,
variants, clones, parental relationships, variant types, genomic positions,
and optional time-point information.
}
\examples{
app <- clevRactApp()

if (interactive()) {
    shiny::runApp(app)
}
}
