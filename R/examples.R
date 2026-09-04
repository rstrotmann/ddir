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


#' Open a DDI assessment report template
#'
#' Crate a new DDI report template as R markdown file in the working directory,
#' and open it. The new report contains examplinib data for reference. Requires
#' RStudio for the editor; otherwise the file is written and the path is
#' returned.
#'
#' @param name Destination file. Defaults to `DDI-assessment.qmd` in the
#'   working directory.
#' @param overwrite Overwrite `name` if it already exists.
#'
#' @returns The destination filename, invisibly.
#' @export
#' @examples
#' \dontrun{
#' ddi_report()
#' }
ddi_report <- function(name = "DDI-assessment.qmd", overwrite = FALSE) {
  template <- ddir_example("DDI-report-examplinib.qmd")

  if (file.exists(name) && !isTRUE(overwrite)) {
    stop(name, " already exists. Set overwrite = TRUE to replace it.")
  }

  if (!file.copy(template, name, overwrite = overwrite)) {
    stop("Could not copy the DDI report template to ", name)
  }

  if (rstudioapi::isAvailable()) {
    rstudioapi::navigateToFile(name)
  } else {
    message("Template written to ", normalizePath(name),
            ". Open it in RStudio or Quarto.")
  }

  invisible(name)
}
