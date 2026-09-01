#' DDI risk object constructor alias
#'
#' @param table The risk table.
#' @param precipitant The precipitant.
#' @param title The title for printing.
#'
#' @returns A risk object.
#' @noRd
new_risk <- function(table, precipitant, title = NULL) {
  structure(
    as_tibble(table),
    class = unique(c("risk", "tbl_df", "tbl", "data.frame")),
    precipitant = precipitant,
    title = title
  )
}


#' DDI risk object constructor alias
#'
#' @param table The risk table.
#' @param precipitant The precipitant.
#' @param title The title for printing.
#'
#' @returns A risk object.
#' @export
risk <- function(table, precipitant, title = NULL) {
  # input validation
  validate_df_argument(table)
  validate_argument(precipitant, "character")
  validate_argument(title, "character", allow_null = TRUE)

  new_risk(table, precipitant, title)
}


#' @importFrom vctrs vec_restore
#' @exportS3Method vctrs::vec_restore
#' @noRd
vec_restore.risk <- function(x, to, ...) {
  new_risk(x, attr(to, "precipitant"), attr(to, "title"))
}


#' reconstruct risk object
#'
#' @param data The risk table.
#' @param template The attributes.
#'
#' @returns A risk object.
#' @exportS3Method dplyr::dplyr_reconstruct
#' @noRd
dplyr_reconstruct.risk <- function(data, template) {
  new_risk(
    data,
    precipitant = attr(template, "precipitant"),
    title = attr(template, "title")
  )
}


#' Print risk object
#'
#' @param x The risk object.
#' @param ... Further arguments.
#'
#' @returns Nothing or a markdodwn-formatted table
#' @exportS3Method base::print
print.risk <- function(x, ...) {
  caption <- attr(x, "title")
  precipitant <-  attr(x, "precipitant")

  x <- mutate(x, across(any_of(
    c("r", "aucr", "Ag", "Ah", "Bg", "Bh", "Cg", "Ch")),
    function(i) round(i, 2)))

  if (isTRUE(getOption("knitr.in.progress"))) {
    caption <- ifelse(
      is.null(caption),
      paste0("DDI risk for precipitant ", precipitant),
      caption
    )
    x |>
      knitr::kable(caption = caption)
  } else {
    cat(paste(hline(), "Clinical DDI risk assessment",  hline(), "\n"))

    if(!is.null(caption))
      cat(paste(caption, "\n\n"))

    cat(df_to_string(x))
  }
}

