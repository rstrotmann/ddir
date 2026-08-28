#' DDI risk object constructor
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

  # business logic
  out <- table
  class(out) <- c("risk", "data.frame")
  attr(out, "precipitant") <- precipitant
  attr(out, "title") <- title

  out
}


#' Print risk object
#'
#' @param x The risk object.
#'
#' @returns Nothing or a markdodwn-formatted table
#' @export
print.risk <- function(x) {
  caption <- attr(x, "title")
  x <- mutate(x, across(any_of("r"), function(i) round(i, 3)))

  if (isTRUE(getOption("knitr.in.progress"))) {
    caption <- ifelse(
      is.null(caption),
      paste0("DDI risk for precipitant ", attr(x, "precipitant")),
      caption
    )
    x$table |>
      knitr::kable(caption = caption)
  } else {
    cat(paste(hline(), "Clinical DDI risk assessment",  hline(), "\n"))

    if(!is.null(caption))
      cat(paste(caption, "\n\n"))

    cat(df_to_string(x))
  }
}

