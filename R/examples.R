#' Get example data path
#'
#' @param path Subdirectory as character.
#'
#' @returns Character.
#' @export
ddir_example <- function(path = NULL) {
  if (is.null(path)) {
    dir(system.file("extdata", package = "ddir"))
  } else {
    system.file("extdata", path, package = "ddir", mustWork = TRUE)
  }
}


#' Open DDI assessment report template
#'
#' @returns Nothing
#' @export
ddi_report <- function() {
  filename <- ddir_example("DDI-report-examplinib.qmd")
  temp <- readChar(filename, file.info(filename)$size)
  rstudioapi::documentNew(temp, type = "rmarkdown")
}
