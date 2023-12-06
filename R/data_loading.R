#' Load perpetrator file
#'
#' @details
#'
#' The file is expected to have the following content:
#'   * Comments start with '#' and go until the end of the line
#'   * lines have the format: __drug__, __parameter__, __value__, __source__
#'   * fields within a line are comma-separated
#'   * fields are plain text without enclosing ""
#'   * the source field is only for reference and usually refers to a clinical
#'     or non-clinical study report. The source field may remain empty. Note
#'     that the line then ends with a comma.
#'
#'   At least the following parameters are expected in the compound file
#'
#'   | parameter | values | detail |
#'   | --------- | ------ | ------ |
#'   | type      | 'parent' or 'metabolite' | the type of substance |
#'   | mw        | numeric | molar weight in g/mol      |
#'   | dose      | numeric | clinical dose in mg        |
#'   | imaxss    | numeric | clinical steady-state Cmax in ng/ml |
#'   | fu        | numeric | unbound fraction in plasma |
#'   | fumic.    | numeric | unbound fraction in microsomes |
#'   | rb        | numeric | blood to plasma ratio      |
#'   | fa        | numeric | fraction absorbed          |
#'   | fg        | numeric | fraction escaping gut metabolism  |
#'   | ka        | numeric | absorption rate constant in 1/min |
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
#'
#' @md
#' @param filename The full path to the compound file
#'
#' @return A list of perpetrator objects
#' @export
#' @import dplyr
#' @import utils
#' @seealso [read_string_perpetrators()]
load_perpetrators <- function(filename) {
  raw <- as.data.frame(read.csv(filename,
                                col.names=c("name", "param", "value", "source"),
                                header = F,
                                comment.char = '#')) %>%
    dplyr::mutate(across(everything(), trimws)) %>%
    dplyr::group_by(name) %>%
    dplyr::group_modify(~ tibble::add_row(param="name", value=.y$name, source="", .x, , .before=1)) %>%
    dplyr::ungroup() %>%
    as.data.frame()

  data <- split(raw, raw$name)
  out <- lapply(data, new_perpetrator)
  return(out)
}


#' Load perpetrator data from csv-formatted string
#'
#' @param input_string The string variable holding the csv-formatted data.
#'
#' @return A list of perpetrator objects
#' @export
#' @import dplyr
#' @import utils
#' @seealso [load_perpetrators()]
read_string_perpetrators <- function(input_string) {
  raw <- as.data.frame(read.csv(text=input_string,
                                col.names=c("name", "param", "value", "source"),
                                header = F,
                                comment.char = '#')) %>%
    dplyr::mutate(across(everything(), trimws)) %>%
    dplyr::group_by(name) %>%
    dplyr::group_modify(~ tibble::add_row(param="name", value=.y$name, source="", .x, , .before=1)) %>%
    dplyr::ungroup() %>%
    as.data.frame()

  data <- split(raw, raw$name)
  out <- lapply(data, new_perpetrator)
  return(out)
}


#' Read perpetrator information from file or text string
#'
#' @param source The connection to read from.
#'
#' @return A list of perpetrator objects.
#' @export
#'
#' @examples
#' read_perpetrators(paste0(system.file(package="ddir"), "/tests/testthat/fixtures/compounds.csv"))
read_perpetrators <- function(source) {
  raw <- as.data.frame(read.csv(source,
                                col.names=c("name", "param", "value", "source"),
                                header = F,
                                comment.char = '#')) %>%
    dplyr::mutate(across(everything(), trimws)) %>%
    dplyr::group_by(name) %>%
    dplyr::group_modify(~ tibble::add_row(param="name", value=.y$name, source="", .x, , .before=1)) %>%
    dplyr::ungroup() %>%
    as.data.frame()

  data <- split(raw, raw$name)
  out <- lapply(data, new_perpetrator)
  return(out)
}


#' Load DMPK data (CYP inhibition, transporter inhibition, etc.)
#'
#' @param filename The filename.
#'
#' @return A list of data frames.
#' @import dplyr
#' @export
load_dmpk_data <- function(filename) {
  raw <- as.data.frame(read.csv(filename,
                                col.names=c("name", "param", "value", "source"),
                                header = F,
                                comment.char = '#')) %>%
    dplyr::mutate(across(everything(), trimws)) %>%
    dplyr::group_by(name) %>%
    dplyr::group_modify(~ tibble::add_row(param="name", value=.y$name, source="", .x, , .before=1)) %>%
    dplyr::ungroup() %>%
    as.data.frame()
  return(raw)
}

#' Load CYP inhibitor data
#'
#' @param filename The filename as string.
#' @return The data as data frame.
#' @export
load_cyp_inhibitor_data <- function(filename) {
  return(load_dmpk_data(filename))
}



#' Read csv-formatted DMPK data from string
#'
#' @param x The string containing the data in csv format.
#' @return The data as data frame.
#' @export
read_string_dmpk_data <- function(x) {
  raw <- as.data.frame(read.csv(text=x,
                                col.names=c("name", "param", "value", "source"),
                                header = F,
                                comment.char = '#')) %>%
    dplyr::mutate(across(everything(), trimws)) %>%
    dplyr::group_by(name) %>%
    dplyr::group_modify(~ tibble::add_row(param="name", value=.y$name, source="", .x, , .before=1)) %>%
    dplyr::ungroup() %>%
    as.data.frame()
  return(raw)
}


#' Read csv-formatted DMPK data
#'
#' @param source The connection to read from.
#'
#' @return The data as data frame.
#' @export
#' @examples
#' read_dmpk(paste0(system.file(package="ddir"), "/tests/testthat/fixtures/cyp-inhibition.csv"))
#'
read_dmpk <- function(source) {
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


#' Load CYP inducer data from file
#'
#' @param filename The filename.
#' @return A list of data frames.
#' @export
load_cyp_inducer_data <- function(filename) {
  raw <- as.data.frame(read.csv(filename,
                                col.names=c("name", "cyp", "emax", "ec50",
                                            "maxc", "source"),
                                header = F,
                                comment.char = '#')) %>%
    dplyr::mutate(across(everything(), trimws)) %>%
    mutate(across(3:5, as.num)) %>%
    as.data.frame()
  return(raw)
}


#' Read csv-formatted CYP inducer data from string
#'
#' @param x The string containing the data in csv format.
#' @return The data as data frame.
#' @export
read_string_cyp_inducer_data <- function(x) {
  raw <- as.data.frame(read.csv(text=x,
                                col.names=c("name", "cyp", "emax", "ec50",
                                            "maxc", "source"),
                                header = F,
                                comment.char = '#')) %>%
    dplyr::mutate(across(everything(), trimws)) %>%
    mutate(across(3:5, as.num)) %>%
    as.data.frame()
  return(raw)
}


#' Read csv-formatted CYP inducer data
#'
#' @param source The connection.
#'
#' @return The inducer data as data frame.
#' @export
#'
#' @examples
#' read_cyp_inducer(paste0(system.file(package="ddir"), "/tests/testthat/fixtures/cyp-induction.csv"))
read_cyp_inducer <- function(source) {
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
