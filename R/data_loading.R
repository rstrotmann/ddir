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
#'   be filled with `NA`. For metabolites, 'oral' must be 'FALSE'.
#'
#'   Note that multiple compounds, e.g., the parent and metabolites may be
#'   included in the perpetrator file. A typical compound file could look e.g.,
#'   like this:
#'
#'   \preformatted{
#'   # name, param, value, source
#'   # parent
#'
#'   examplinib,  oral,     TRUE,
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
#'   M1,  oral,   FALSE,
#'   M1,  mw,     506.56,
#'   M1,  dose,   NA,
#'   M1,  imaxss, 1038,      study 001
#'   M1,  fu,     0.012,     study 002
#'   M1,  fumic,  1,         default
#'   M1,  rb,     1,         study 002
#'   M1,  fa,     NA,
#'   M1,  fg,     NA,
#'   M1,  ka,     NA,
#'   }
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


#' Read csv-formatted CYP inhibition data
#'
#' This function loads CYP inhibition data from a csv file. The expected fields
#' are (in this order) the compound name, the CYP enzyme, the Ki and the source
#' information for the respective data. The latter field may remain empty.
#'
#' Comment lines must start with '#'.
#'
#' @details
#' A valid source is, e.g.,
#' \preformatted{
#' # PARENT
#' # name,     CYP,     Ki,   source
#' examplinib, CYP1A2,  NA,
#' examplinib, CYP2B6,  NA,
#' examplinib, CYP2C8,  11,   study 001
#' examplinib, CYP2C9,  13.5, study 001
#' examplinib, CYP2C19, 15,   study 001
#' examplinib, CYP2D6,  NA,
#' examplinib, CYP3A4,  12.5, study 001

#' # METABOLITE
#' # name,     CYP,     Ki,   source
#' M1,         CYP2C9,  4.4,  study 002
#' }
#'
#' @param source The connection to read from.
#'
#' @return A data frame.
#' @export
#' @examples
#' read_inhibitor_data(textConnection(examplinib_cyp_inhibition_string))
read_inhibitor_data <- function(source) {
  raw <- as.data.frame(read.csv(source,
                                col.names=c("name", "item", "ki", "source"),
                                header = F,
                                blank.lines.skip = TRUE,
                                comment.char = '#')) %>%
    dplyr::mutate(across(everything(), trimws)) %>%
    dplyr::filter(name != "") %>%
    as.data.frame()
  return(raw)
}


#' Read UGT inhibition data
#'
#' Read UGT inhibition data from a file or text connection.
#' @details
#' The following, comma-separated fields are expected in the input:
#' * 'name' The perpetrator compound name
#' * 'ugt' The UGT enzyme as (upper case) character.
#' * 'ic50' The \eqn{IC_{50}}
#' * 'source' Optional source information as character.
#'
#' Lines starting with '#' are considered comments and are not evaluated.
#' @details
#' The following is an example of a valid input:
#' \preformatted{
#' # PARENT
#' # compound, enzyme, IC50, source
#' examplinib, UGT1A1, 15, study 009
#' examplinib, UGT1A3, 15, study 009
#' examplinib, UGT1A4, 15, study 009
#' examplinib, UGT1A6, 15, study 009
#' examplinib, UGT1A9, 3.8, study 009
#' examplinib, UGT2B7, 15, study 009
#' examplinib, UGT2B15, 15, study 009
#' examplinib, UGT2B17, 6.1, study 009
#' # METABOLITE
#' # compound, enzyme, IC50, source
#' M1, UGT1A1, 1.1, study 009
#' M1, UGT1A3, 5.8, study 009
#' M1, UGT1A4, 6.2, study 009
#' M1, UGT1A6, 15, study 009
#' M1, UGT1A9, 3.6, study 009
#' M1, UGT2B7, 15, study 009
#' M1, UGT2B15, 9.6, study 009
#' }
#' @param source The file or text connection to read from.
#'
#' @return A data frame.
#' @export
#'
#' @examples
#' read_ugt_inhibitor_data(textConnection(examplinib_ugt_inhibition_string))
read_ugt_inhibitor_data <- function(source) {
  raw <- read_inhibitor_data(source)
  colnames(raw) <- c("name", "ugt", "ic50", "source")
  return(raw)
}


#' Read CYP inhibition data
#'
#' Read CYP inhibition data from a file or text connection.
#' @details
#' The following, comma-separated fields are expected (in this order):
#' * 'name' The perpetrator compound name
#' * 'cyp' The UGT enzyme as (upper case) character.
#' * 'ki' The \eqn{k_i}
#' * 'source' Optional source information as character.
#'
#' Lines starting with '#' are considered comments and are not evaluated.
#' @details
#' The following is an example of a valid input:
#' \preformatted{
#' # PARENT
#' # compound, CYP, ki, source
#' examplinib, CYP1A2,  NA,
#' examplinib, CYP2B6,  NA,
#' examplinib, CYP2C8,  11,   study 001
#' examplinib, CYP2C9,  13.5, study 001
#' examplinib, CYP2C19, 15,   study 001
#' examplinib, CYP2D6,  NA,
#' examplinib, CYP3A4,  12.5, study 001
#' # METABOLITE
#' M1,         CYP2C9,  4.4,  study 002
#' }
#' @param source The file or text connection to read from.
#' @return A data frame.
#' @export
#' @noRd
#' @examples
#' read_cyp_inhibitor_data(textConnection(examplinib_cyp_inhibition_string))
read_cyp_inhibitor_data <- function(source) {
  raw <- read_inhibitor_data(source)
  colnames(raw) <- c("name", "cyp", "ki", "source")
  return(raw)
}


#' Read csv-formatted CYP TDI data
#'
#' This function reads comma-separated CYP TDI data from a file or text
#' connection.
#' @details
#' The following fields are expected (in this order):
#' * 'name' The perpetrator compound name as character.
#' * 'cyp' The CYP enzyme as character.
#' * 'ki' The \eqn{K_I} in \eqn{\mu M} as numeric.
#' * 'kinact' The \eqn{k_{inact}} in \eqn{1/h} as numeric.
#' * 'source' Optional source information as character.
#' Lines starting with '#' are interpreted as comments and are not evaluated.
#'
#' A valid data set is, e.g.,
#' \preformatted{
#' # compound, CYP,    ki,   kinact, source
#' examplinib, CYP3A4, 0.17, 0.04, study 001
#' }
#'
#' @param source The connection to read from.
#'
#' @return A data frame.
#' @export
read_tdi_data <- function(source) {
  raw <- as.data.frame(read.csv(source,
                                col.names=c("name", "cyp", "ki", "kinact", "source"),
                                header = F,
                                blank.lines.skip = TRUE,
                                comment.char = '#')) %>%
    dplyr::mutate(across(everything(), trimws)) %>%
    dplyr::filter(name != "") %>%
    mutate(across(3:4, as.num)) %>%
    as.data.frame()
  return(raw)
}


#' Read csv-formatted CYP inducer data
#'
#' This function loads CYP inducer data from a csv file. The expected fields are
#' (in this order) the compound name, the CYP enzyme, the Emax, the EC50, the
#' maximal tested concentration and the source reference. The latter field may
#' remain empty.
#'
#' Comment lines must start with '#'.
#'
#' @details
#' A valid source is, e.g.,
#' \preformatted{
#' # PARENT
#' # compound, CYP, Emax, EC50, max c, source
#' examplinib, CYP1A2,  1,    NA,   5,  study 007
#' examplinib, CYP2B6,  1,    NA,   5,  study 007
#' examplinib, CYP2C8,  NA,   NA,   NA,
#' examplinib, CYP2C9,  NA,   NA,   NA,
#' examplinib, CYP2C19, NA,   NA,   NA,
#' examplinib, CYP2D6,  NA,   NA,   NA,
#' examplinib, CYP3A4,   7.35, 1.64, 3,  study 007
#'
#' # METABOLITE
#' # compound, CYP, ki, source
#' M1, CYP1A2,  1,    NA,   5,  study 007
#' M1, CYP2B6,  6.98, 1.86, 5,  study 007
#' M1, CYP2C8,  NA,   NA,   NA,
#' M1, CYP2C9,  NA,   NA,   NA,
#' M1, CYP2C19, NA,   NA,   NA,
#' M1, CYP2D6,  NA,   NA,   NA,
#' M1, CYP3A4,  22.7, 1.1,  5,  study 007
#' }
#'
#' @param source The connection to read from.
#' @return The data as data frame.
#' @export
#' @examples
#' read_inducer_data(textConnection(examplinib_cyp_induction_string))
#'
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


#' Read transporter inhibition data
#'
#' Read transporter inhibition data from a file or text connection.
#' @details
#' The following, comma-separated fields are expected (in this order):
#' * 'name' The perpetrator compound name
#' * 'cyp' The UGT enzyme as (upper case) character.
#' * 'ki' The \eqn{k_i}
#' * 'source' Optional source information as character.
#'
#' Lines starting with '#' are considered comments and are not evaluated.
#' @details
#' The following is an example of a valid input:
#' \preformatted{
#' # PARENT
#' # name,     transporter, IC50, source
#' examplinib, Pgp,         0.41, study 005
#' examplinib, BCRP,        1.9,  study 005
#' examplinib, OCT1,        2.3,  study 006
#' examplinib, OATP1B1,     177,  study 006
#' examplinib, OATP1B3,     35,   study 006
#' examplinib, OAT1,        271,
#' examplinib, OAT3,        300,
#' examplinib, BSEP,        12.8,
#' examplinib, OCT2,        67,   study 006
#' examplinib, MATE1,       3.6,  study 006
#' examplinib, MATE2k,      1.1,  study 006
#' }
#' @param source The file or text connection to read from.
#' @return A data frame.
#' @export
#'
#' @examples
#' read_transporter_inhibitor_data(textConnection(examplinib_transporter_inhibition_string))
read_transporter_inhibitor_data <- function(source) {
  raw <- read_inhibitor_data(source)
  colnames(raw) <- c("name", "transporter", "ic50", "source")
  return(raw)
}

