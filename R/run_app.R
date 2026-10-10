#' Launch the iOmicsPASS+ Shiny GUI
#'
#' Opens the interactive graphical user interface for iOmicsPASS+ in a web
#' browser.  The GUI guides users through data upload, optional network
#' inference, subnetwork signature analysis, and prediction on new samples.
#'
#' @param launch.browser Logical.  If \code{TRUE} (default) the system default
#'   browser is opened automatically.  Set to \code{FALSE} to start the server
#'   without opening a browser (useful when running on a remote server).
#' @return Called for its side effect; does not return a value.
#' @examples
#' \dontrun{
#' run_iOmicsPASS()
#' }
#' @export
run_iOmicsPASS <- function(launch.browser = TRUE) {
  app_dir <- system.file("shiny", "iOmicsPASSplus", package = "iOmicsPASSplus")
  if (app_dir == "") {
    stop("Shiny app directory not found. Please re-install iOmicsPASSplus.",
         call. = FALSE)
  }
  shiny::runApp(app_dir, launch.browser = launch.browser)
}
