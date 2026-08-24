#' Generic DDI risk class definition
#'
#' @slot precipitant Name of perpetrator.
#' @slot table Tabulated risk assessment data.
#'
#' @returns An ddi_risk class object.
#' @export
setClass(
  Class = "ddi_risk",
  representation(
    precipitant = "perpetrator",
    table = "data.frame",
    title = "character"
  ),
  prototype(
    precipitant = new("perpetrator"),
    table = data.frame(),
    title = ""
  )
)


#' Show ddi_risk object
#'
#' @param ddi_risk Ddi_risk object.
#' @returns Nothing.
#' @export
setMethod(
  "show", "ddi_risk",
  function(object) {
    caption <- ifelse(
      object@title == "",
      paste0("DDI risk for precipitant ", object@precipitant@name),
      object@title
    )
    temp <- object@table
    cat(paste0(caption, "\n\n"))
    cat(df_to_string(object@table))
  }
)


#' Print ddi_risk object
#'
#' @param ddi_risk Ddi_risk object.
#' @returns Markdown-formatted text.
#' @export
setMethod(
  "print", "ddi_risk",
  function(x) {
    temp <- x@table
    caption <- ifelse(
      x@title == "",
      paste0("DDI risk for precipitant ", x@precipitant@name),
      x@title
    )

    x@table |>
      mutate(
        across(any_of(c(
          "r", "r_gut", "aucr", "Ag", "Ah", "Bg", "Bh", "Cg", "Ch")),
        function(r) round(r, 2))) |>
      mutate(across(
        starts_with("risk"),
        function(x) case_when(
          is.na(x) ~ "",
          x == TRUE ~ "Yes",
          .default = "No"
        )
      )) |>
      knitr::kable(caption = caption, col.names = make_labels(x@table))
  }
)
