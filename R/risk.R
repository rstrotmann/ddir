#' Generic DDI risk class definition
#'
#' @slot object Name of perpetrator.
#' @slot table Tabulated risk assessment data.
#'
#' @returns An ddi_risk class object.
#' @export
setClass(
  Class = "ddi_risk",
  representation(
    object = "perpetrator",
    table = "data.frame",
    title = "character"
  ),
  prototype(
    object = new("perpetrator"),
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
      paste0("DDI risk for object ", object@object@name),
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
      paste0("DDI risk for object ", x@object@name),
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
      kable(caption = caption, col.names = make_labels(x@table))
  }
)
