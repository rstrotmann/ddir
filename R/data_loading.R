#' Read perpetrator information from file or text string
#'
#' The input source can either be a file name as string or a text connection.
#' This can be used to read the compound data from a string.
#'
#' @details
#' The following fields are expected in the input (in this order):
#' * 'name' The compound name as character.
#' * 'parameter' The parameter as character.
#' * 'value' The respective value as character.
#' * 'source' Optional source information as character.
#'
#' The following parameters are expected:
#' * 'oral' A Boolean (i.e., `TRUE` or `FALSE`) to indicate whether the drug is
#' subject to first-pass effects.
#' * 'mw' The molar weight in g/mol as numeric.
#' * 'dose' The clinically administered dose in mg as numeric.
#' * 'imaxss' The (total) \eqn{C_{max}} in ng/ml after administration of the
#' clinical dose.
#' * 'fu' The free unbound) fraction of the drug in plasma.
#' * 'fumic' The free (unbound) fraction in microsomal preparations.
#' * 'rb' The blood-to-plasma concentration ratio.
#' * 'fa' The fraction absorbed of the drug.
#' * 'fg' The fraction of the administered dose escaping gut metabolism.
#' * 'ka' The absorption rate constant in 1/min.
#' * 'solubility' The aqueous solubility of the compound in mg/l.
#'
#' Lines starting with '#' are considered comments and are not evaluated.
#'
#' Note that multiple compounds, e.g., the parent and metabolites may be
#' included in the perpetrator file. The below is an example of a valid compound
#' file:
#'
#' \preformatted{
#' # name, param, value, source
#' # parent
#' examplinib,  oral,     TRUE,
#' examplinib,  mw,       492.6,
#' examplinib,  dose,     450,       clinical dose
#' examplinib,  imaxss,   3530,      study 001
#' examplinib,  fu,       0.023,     study 002
#' examplinib,  fumic,    1,         default
#' examplinib,  rb,       1,         study 003
#' examplinib,  fa,       0.81,      study 003
#' examplinib,  fg,       1,         default
#' examplinib,  ka,       0.00267,   unknown
#'
#' # metabolite
#' M1,  oral,   FALSE,
#' M1,  mw,     506.56,
#' M1,  dose,   NA,
#' M1,  imaxss, 1038,      study 001
#' M1,  fu,     0.012,     study 002
#' M1,  fumic,  1,         default
#' M1,  rb,     1,         study 002
#' M1,  fa,     NA,
#' M1,  fg,     NA,
#' M1,  ka,     NA,
#' }
#'
#' @param source The file name or text connection to read from.
#' @return A perpetrator object if only one compound in the input source, or
#' list of perpetrator objects.
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
  if(length(out) == 1) {
    return(out[[1]])
  } else {
    return(out)
  }
}


#' Read csv-formatted CYP inhibition data
#'
#' This function loads CYP inhibition data from a csv file. The expected fields
#' are (in this order) the compound name, the CYP enzyme, the Ki and the source
#' information for the respective data. The latter field may remain empty.
#'
#' Comment lines must start with '#'.
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is deprecated in favor of [read_cyp_inhibitor_data()],
#' [read_ugt_inhibitor_data()] or [read_transporter_inhibitor_data()].
#' @details
#' A valid source is, e.g.,
#'
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
#' * 'name' The perpetrator compound name as character.
#' * 'ugt' The UGT enzyme as (upper case) character.
#' * 'ic50' The \eqn{IC_{50}} in µM as numeric.
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
#' The following, comma-separated fields are expected in the input (in this
#' order):
#' * 'name' The perpetrator compound name as character.
#' * 'cyp' The UGT enzyme as (upper case) character.
#' * 'ki' The \eqn{k_i} in µM as numeric.
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
#' * 'cyp' The CYP enzyme as (upper case) character.
#' * 'ki' The \eqn{K_I} in µM as numeric.
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
#' This function loads CYP inducer data from a csv file. The following fields
#' are expected:
#' * 'name' The name of the perpetrator compound as character.
#' * 'cyp' The CYP enzyme as (upper case) character.
#' * 'emax' The \eqn{E_{max}}, i.e., the maximum induction effect determined in
#' vitro as numeric.
#' * 'ec50' The \eqn{EC_{50}} in µM as numeric.
#' * 'maxc' The maximal concentration in µM tested in the in vitro assay as
#' numeric.
#' * 'source' Optional source information as character.
#'
#' Comment lines must start with '#'.
#' @details
#' A valid source is, e.g.,
#' \preformatted{
#' # PARENT
#' # compound, CYP,     Emax, EC50, max c, source
#' examplinib, CYP1A2,  1,    NA,   5,     study 007
#' examplinib, CYP2B6,  1,    NA,   5,     study 007
#' examplinib, CYP2C8,  NA,   NA,   NA,
#' examplinib, CYP2C9,  NA,   NA,   NA,
#' examplinib, CYP2C19, NA,   NA,   NA,
#' examplinib, CYP2D6,  NA,   NA,   NA,
#' examplinib, CYP3A4,   7.35, 1.64, 3,    study 007
#'
#' # METABOLITE
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
#' * 'name' The perpetrator compound name as character.
#' * 'cyp' The UGT enzyme as (upper case) character.
#' * 'ic50' The \eqn{IC_{50}} of the inhibition in μM.
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
#' @examples
#' read_transporter_inhibitor_data(textConnection(examplinib_transporter_inhibition_string))
read_transporter_inhibitor_data <- function(source) {
  raw <- read_inhibitor_data(source)
  colnames(raw) <- c("name", "transporter", "ic50", "source")
  return(raw)
}

