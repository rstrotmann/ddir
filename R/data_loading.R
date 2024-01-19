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
#' ```{r echo=F, comment=NA}
#' cat(examplinib_compounds_string)
#' ```
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
    filter(trimws(name) != "") %>%
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
#' ```{r echo=F, comment=NA}
#' cat(examplinib_cyp_inhibition_string)
#' ```
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
#' ```{r echo=F, comment=NA}
#' cat(examplinib_ugt_inhibition_string)
#' ```
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
#' * 'cyp' The CYP enzyme as (upper case) character.
#' * 'ki' The \eqn{k_i} in µM as numeric.
#' * 'source' Optional source information as character.
#'
#' Lines starting with '#' are considered comments and are not evaluated.
#' @details
#' The following is an example of a valid input:
#' ```{r echo=F, comment=NA}
#' cat(examplinib_cyp_inhibition_string)
#' ```
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
#' ```{r echo=F, comment=NA}
#' cat(examplinib_cyp_tdi_string)
#' ```
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
#' ```{r echo=F, comment=NA}
#' cat(examplinib_cyp_induction_string)
#' ```
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
#' ```{r echo=F, comment=NA}
#' cat(examplinib_transporter_inhibition_string)
#' ```
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

