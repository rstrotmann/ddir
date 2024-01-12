#' Read perpetrator information from file or text string
#'
#' The input source can either be a file name as string or a text connection.
#' This can be used to read the compound data from a string.
#'
#' @details
#'
#' The file is expected to comply with the following:
#'
#'   * Comments start with '#' and go until the end of the line
#'
#'   * lines have the format: __drug__, __parameter__, __value__, __source__
#'
#'   * fields within a line are comma-separated
#'
#'   * fields are plain text without enclosing ""
#'
#'   * the source field is only for reference and usually refers to a clinical
#'     or non-clinical study report. The source field may remain empty. Note
#'     that the line then ends with a comma.
#'
#'
#'   For metabolites, the fields dose, fa, fg and ka are not relevant and should
#'   be filled with `NA`.
#'
#'   Note that multiple compounds, e.g., the parent and metabolites may be
#'   included in the perpetrator file. A typical compound file could look e.g.,
#'   like this:
#'
#'   ```
#'   # name, param, value, source
#'   # parent
#'
#'   examplinib,  type,     parent,
#'   examplinib,  mw,       492.6,
#'   examplinib,  dose,     450,       clinical dose
#'   examplinib,  imaxss,   3530,      study 001
#'   examplinib,  fu,       0.023,     study 002
#'   examplinib,  fumic,    1,         default
#'   examplinib,  rb,       1,         study 003
#'   examplinib,  fa,       0.81,      study 003
#'   examplinib,  fg,       1,         default
#'   examplinib,  ka,       0.00267,   unknown
#'
#'   # metabolite
#'
#'   M1,  type,   metabolite,
#'   M1,  mw,     506.56,
#'   M1,  dose,   NA,
#'   M1,  imaxss, 1038,      study 001
#'   M1,  fu,     0.012,     study 002
#'   M1,  fumic,  1,         default
#'   M1,  rb,     1,         study 002
#'   M1,  fa,     NA,
#'   M1,  fg,     NA,
#'   M1,  ka,     NA,
#'   ```
#' @param source The connection to read from.
#'
#' @return A list of perpetrator objects.
#' @import dplyr
#' @export
#'
#' @examples
#' read_perpetrators(textConnection(examplinib_compounds_string))
read_perpetrators <- function(source) {
  raw <- as.data.frame(read.csv(source,
                                col.names=c("name", "param", "value", "source"),
                                header = F,
                                comment.char = '#')) %>%
    dplyr::mutate(across(everything(), trimws)) %>%
    dplyr::group_by(name) %>%
    dplyr::group_modify(~ tibble::add_row(param="name", value=.y$name,
                                          source="", .x, , .before=1)) %>%
    dplyr::ungroup() %>%
    as.data.frame()

  data <- split(raw, raw$name)
  out <- lapply(data, function(x) {x <- new_perpetrator(x %>% select(-name))})
  return(out)
}


#' Read csv-formatted DMPK data
#'
#' @param source The connection to read from.
#'
#' @return The data as data frame.
#' @export
read_inhibitor_data <- function(source) {
  raw <- as.data.frame(read.csv(source,
                                col.names=c("name", "param", "value", "source"),
                                header = F,
                                comment.char = '#')) %>%
    dplyr::mutate(across(everything(), trimws)) %>%
    dplyr::group_by(name) %>%
    dplyr::group_modify(~ tibble::add_row(param="name",
                                          value=.y$name,
                                          source="",
                                          .x, , .before=1)) %>%
    dplyr::ungroup() %>%
    as.data.frame()
  return(raw)
}


#' #' Read CYP inhibitor data
#' #'
#' #' @param source The connection to read from.
#' #'
#' #' @return The data as data frame.
#' #' @export
#' read_cyp_inhibitor_data <- function(source) {
#'   return(read_dmpk_data(source))
#' }


#' Read csv-formatted CYP inducer data
#'
#' @param source The connection to read from.
#' @return The data as data frame.
#' @export
read_inducer_data <- function(source) {
  raw <- as.data.frame(read.csv(source,
                                col.names=c("name", "cyp", "emax", "ec50",
                                            "maxc", "source"),
                                header = F,
                                comment.char = '#')) %>%
    dplyr::mutate(across(everything(), trimws)) %>%
    mutate(across(3:5, as.num)) %>%
    as.data.frame()
  return(raw)
}



