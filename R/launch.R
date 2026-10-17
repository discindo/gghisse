#' Launch the gghisse Shiny web application
#'
#' @return None, called for the side effect of running the 'Shiny'
#'   application.
#' @examples
#' if (interactive()) {
#'   launch_gghisse_app()
#' }
#' @export
launch_gghisse_app <- function() {
  app_dir <- system.file("gghisse-web", "gghisse", package = "gghisse")
  if (app_dir == "") {
    stop(
      "Could not find example directory. Try re-installing `gghisse`.",
      call. = FALSE
    )
  }

  shiny::runApp(app_dir, display.mode = "normal")
}
