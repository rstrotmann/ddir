
convert_df_to_perp <- function(x) {
  out <- perpetrator("", TRUE, 0, 0, 0)
  for (p in c("mw", "dose", "imaxss", "fu", "fummic", "rb", "fa", "fg", "ka", "solubility")) {
    if (p %in% x$param) {
      slot(out, p) <- as.numeric(x[x$param == p, "value"])
    }
  }
  if ("name" %in% x$param) slot(out, "name") <- x[x$param == "name", "value"]
  if ("oral" %in% x$param) slot(out, "oral") <- as.logical(x[x$param == "oral", "value"])
  source <- x$source
  names(source) <- x$param
  methods::slot(out, "source") <- source[source != ""]
  out
}





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
#' @param source The file name or text connection to read from.
#' @return A perpetrator object if only one compound in the input source, or
#' list of perpetrator objects.
#' @import dplyr
#' @noRd
read_perpetrators <- function(source) {
  raw <- as.data.frame(read.csv(
    source,
    col.names=c("name", "param", "value", "source"),
    header = F,
    comment.char = '#')
  ) |>
    filter(trimws(name) != "") |>
    dplyr::mutate(across(everything(), trimws)) |>
    dplyr::group_by(name) |>
    dplyr::group_modify(~ tibble::add_row(param="name", value=.y$name,
                                          source="", .x, , .before=1)) |>
    dplyr::ungroup() |>
    as.data.frame()

  data <- split(raw, raw$name)
  out <- lapply(
    data,
    function(x) {
      x <- convert_df_to_perp(x)
    }
  )

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
#' @noRd
read_inhibitor_data <- function(source) {
  raw <- as.data.frame(read.csv(
    source,
    col.names=c("name", "item", "ki", "source"),
    header = F,
    blank.lines.skip = TRUE,
    comment.char = '#')
  ) %>%
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
#' @noRd
read_ugt_inhibitor_data <- function(source) {
  raw <- read_inhibitor_data(source)
  colnames(raw) <- c("name", "object", "ic50", "source")

  data <- split(raw, raw$name)
  out <- lapply(
    data,
    function(x) {
      x <- x |>
        # rename(object = ugt) |>
        mutate(ic50 = as.numeric(ic50)) |>
        filter(!is.na(ic50))
      inhibitor(data = select(x, -name), object <- unique(x$name))
    }
  )

  if(length(out) == 1) {
    return(out[[1]])
  } else {
    return(out)
  }
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
#' @noRd
read_cyp_inhibitor_data <- function(source) {
  raw <- read_inhibitor_data(source)
  colnames(raw) <- c("name", "object", "ki", "source")
  data <- split(raw, raw$name)
  out <- lapply(
    data,
    function(x) {
      suppressWarnings(
        x <- x |>
          mutate(ki = as.numeric(ki))
      )
        # filter(!is.na(ki))
      inhibitor(data = select(x, -name), precipitant <- unique(x$name))
    }
  )

  if(length(out) == 1) {
    return(out[[1]])
  } else {
    return(out)
  }
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
#' @noRd
read_tdi_data <- function(source) {
  raw <- as.data.frame(read.csv(
    source,
    col.names=c("name", "object", "ki", "kinact", "source"),
    header = F,
    blank.lines.skip = TRUE,
    comment.char = '#')
  ) |>
    dplyr::mutate(across(everything(), trimws)) |>
    dplyr::filter(name != "") |>
    mutate(across(3:4, as.num)) |>
    as.data.frame()

  data <- split(raw, raw$name)
  out <- lapply(
    data,
    function(x) {
      x <- x |>
        mutate(ki = as.numeric(ki)) |>
        filter(!is.na(ki))
      inhibitor(data = select(x, -name), object <- unique(x$name))
    }
  )

  if(length(out) == 1) {
    return(out[[1]])
  } else {
    return(out)
  }

  # out <- inhibitor(select(raw, -name))
  # return(raw)
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
#' @noRd
read_inducer_data <- function(source) {
  raw <- as.data.frame(read.csv(
    source,
    col.names=c("name", "cyp", "emax", "ec50",
                "maxc", "source"),
    header = F,
    comment.char = '#')
  ) |>
    dplyr::mutate(across(everything(), trimws)) |>
    mutate(across(3:5, as.num)) |>
    as.data.frame()

  data <- split(raw, raw$name)
  out <- lapply(
    data,
    function(x) {
      x <- x |>
        mutate(across(c("emax", "ec50", "maxc"), as.numeric)) |>
        rename(object = cyp) |>
        mutate(max_c = maxc) |>
        filter(!is.na(emax))
      inducer(data = select(x, -name), object <- unique(x$name))
    }
  )

  if(length(out) == 1) {
    return(out[[1]])
  } else {
    return(out)
  }
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
#' @noRd
read_transporter_inhibitor_data <- function(source) {
  raw <- read_inhibitor_data(source)
  colnames(raw) <- c("name", "object", "ic50", "source")

  data <- split(raw, raw$name)
  out <- lapply(
    data,
    function(x) {
      x <- x |>
        mutate(ic50 = as.numeric(ic50)) |>
        filter(!is.na(ic50))
      inhibitor(data = select(x, -name), object <- unique(x$name))
    }
  )

  if(length(out) == 1) {
    return(out[[1]])
  } else {
    return(out)
  }
}
