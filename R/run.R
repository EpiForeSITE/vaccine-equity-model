
#' Runs the shiny app
#' @param ... Further arguments passed to [shiny::shinyAppDir()].
#' @details
#' The app featured in this package is the one presented in the shiny demo:
#' <https://shiny.posit.co/r/getstarted/shiny-basics/lesson1/index.html>.
#'
#' @return
#' Starts the execution of the app, printing the port
#' on the console.
#' @export
#' @examples
#' # To be executed interactively only
#' if (interactive()) {
#'   run_my_app()
#' }
#' @import shiny
#' @importFrom htmltools tags img
#' @importFrom bslib page_sidebar sidebar
#' @importFrom graphics hist
#' @importFrom utils data
#' @export
#'
#'
run_my_app <- function(...) {
  shiny::shinyAppDir(
    system.file("app/", package = "vaccine.equity")
  )
}
