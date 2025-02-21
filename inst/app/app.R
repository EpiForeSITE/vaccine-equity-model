#' Runs the vaccine model
#'
#' not needed?


#' @export
run_model <- function() {
  shiny::shinyAppDir(
    system.file("app/", package = "vaccine.equity")
  )
}
