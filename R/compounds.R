#' Render data frame to character.
#'
#' This function renders a data frame into a string object.
#'
#' @param df The data frame.
#' @param indent A string that defines the left indentation of the rendered
#' output.
#' @param colnames Boolean value to indicate whether column names are to be
#' included in the output.
#' @param n The number of lines to be rendered. If NULL (default), all lines
#' are rendered.
#'
#' @return The data frame representation as character.
#' @import stringr
#' @import utils
#' @noRd
df_to_string <- function(df, indent="", n=NULL, colnames=TRUE){
  df <- as.data.frame(df)
  max.widths <- as.numeric(
    lapply(rbind(df, names(df)),
           FUN=function(x) max(sapply(as.character(x), nchar), na.rm=TRUE)))
  line = df[1,]

  render.line <- function(line){
    out <- indent
    for(i in 1:length(line)){
      out <- paste0(out, sprintf(paste0("%-", max.widths[i]+3, "s"),
                                 as.character(line[i])))
    }
    return(out)
  }
  out <- NULL
  if(colnames){
    out <- render.line(data.frame(as.list(names(df))))
  }
  if(!is.null(n)){
    df <- head(df, n=n)
  }
  for(i in 1:nrow(df)){
    out <- paste(out, render.line(df[i,]), sep="\n")
  }
  return(stringr::str_trim(out))
}


#' Convert field to numeric with NA translated to 0
#'
#' @param x The input as character.
#' @param na.strings Strings representing NA values.
#' @return Numeric.
#' @noRd
as.num = function(x, na.strings = "NA") {
  stopifnot(is.character(x))
  na = x %in% na.strings
  x[na] = "0"
  x = as.numeric(x)
  x[na] = NA_real_
  x
}


#' Perpetrator names as comma-separated string
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' @param perps perpetrator objects as a list.
#' @import lifecycle
#' @return The output as string.
#' @export
#' @noRd
names_string <- function(perps) {
  lifecycle::deprecate_warn("0.8.1", "names_string()", "compound_names_string()")
  return(paste(lapply(perps, function(p) {p[(p$param=="name"), "value"]}), collapse=", "))
}


#' Compound names as single string
#'
#' This function returns a single string with the names of the perpetrator
#' compounds listed nicely.
#' @param compounds A list of perpetrator objects
#' @return A character string.
#' @export
#' @examples
#' compound_names_string(examplinib_compounds)
compound_names_string <- function(compounds) {
  temp <- lapply(compounds, function(p) {
    p[(p$param=="name"), "value"]})

  if(length(temp) == 1) {
    return(temp[[1]])
  }

  if(length(temp) > 1) {
    return(paste(paste(temp[1:length(temp)-1], collapse = ", "), "and",
           temp[length(temp)]))
  }
}


#' Perpetrator object constructor
#'
#' This function converts a data frame into a perpetrator object.
#'
#' @details
#' The input is a data frame with the columns 'param', 'value' and source'.
#'
#' Rows with the following parameters ('param') are expected in the input:
#'
#' | PARAM | REQUIREMENT | DESCRIPTION |
#' | --- | --- | --- |
#' | name | required |  |
#' | type | optional | either 'parent' (default) or 'metabolite' |
#' | mw | required | molar weight in g/mol |
#' | dose | required for parent compounds | clinical dose in mg |
#' | imaxss | required | total Cmax in ng/ml |
#' | fu | optional | fraction unbound, default is 1 |
#' | fumic | optional | microsomal unbound fraction, dafault is 1 |
#' | rb | optional | blood-to-plasma concentration ratio, default is 1 |
#' | fa | optional | fraction absorbed, default is 1 |
#' | fg | optional | fraction escaping gut metabolism, default is 1 |
#' | ka | optional | absorption rate constant, default is 0.1 /min |
#' | solubility | optional | solubility of the compound in mg/l, defaults to `Inf` |
#'
#' The following example is an example for a valid input data frame:
#'
#' \preformatted{
#'         param      value        source
#' 1        name examplinib
#' 2        type     parent
#' 3          mw      492.6
#' 4        dose        450 clinical dose
#' 5      imaxss       3530     study 001
#' 6          fu      0.023     study 002
#' 7       fumic          1       default
#' 8          rb          1     study 003
#' 9          fa       0.81     study 003
#' 10         fg          1       default
#' 11         ka    0.00267       unknown
#' 12 solubility        Inf       default
#' }
#'
#' @param df A data frame to be converted into a perpetrator object. See
#' 'Details' for the expected format and fields.
#' @return A perpetrator object.
#' @import assertr
#' @import tibble
#' @import tidyr
#' @seealso [read_perpetrators()]
#' @export
new_perpetrator <- function(df) {
  default_values <- tribble(
    ~param, ~default,
    "name",   NA,
    # "type",   "parent",
    "oral",   "TRUE",
    "mw",     NA,
    "dose",   NA,
    "imaxss", NA,
    "fu",     "1",
    "fumic",  "1",
    "rb",     "1",
    "fa",     "1",
    "fg",     "1",
    "ka",     "0.1",
    "solubility", "Inf"
  )

  out <- df %>%
    assertr::verify(assertr::has_all_names("param", "value", "source")) %>%
    assertr::verify(c("name", "mw", "dose", "imaxss") %in% .$param) %>%
    as.data.frame()

  if("type" %in% (out %>% pull(param))) {
    message(paste("Deprecation warning: The input contains the field 'type'.",
      "This field was expected\nto have the values 'parent' or 'metabolite'",
      "to indicate whether the compound is\nsubject to first-pass metabolism.",
      "This field has been deprecated. Please use instead\n'oral' with the",
      "possible, values of 'TRUE' or 'FALSE'!\nFor compatibility reasons,",
      "the field 'oral' bas been autogenerated and your input\nshould still",
      "work. However, this function may be deleted in future releases."))
    if(!("oral" %in% (out %>% pull(param)))) {
      out <- out %>% add_row(param="oral",
        value=as.character(as.character(df[which(df$param=="type"), "value"])=="parent"),
        source="autogenerated from 'type'")
    }
  }

  out <- out %>%
    right_join(default_values, by="param") %>%
    mutate(source = case_when(is.na(value) & !is.na(default) ~ "default",
                              .default = source)) %>%
    mutate(value = case_when(is.na(value) & !is.na(default) ~ default,
                             .default = value)) %>%
    select(-default) %>%
    assertr::verify(!is.na(value)) %>%
    remove_rownames() %>%
    magrittr::set_rownames(.$param)

  class(out) <- c("perpetrator", "data.frame")
  return(out)
}


#' Make a perpetrator object from compound data
#'
#' This function creates a perpetrator object from key compound data. Note that
#' the 'source' field is empty. If you want to include source information, you
#' can either use [new_perpetrator()] which takes as input a data that may
#' include source information, or you can create a perpetrator object from a
#' text input string using [read_perpetrators()].
#'
#' @param name The compound name as character.
#' @param dose The clinical dose in mg.
#' @param imaxss The (total) steady-state Cmax in ng/ml.
#' @param mw The molar weight in g/mol.
#' @param type The compound type as character. Must be either 'parent' or
#' metabolite'.
#' @param fu The fraction unbound as numeric. Defaults to 1.
#' @param fumic The fraction unbound in the microsomes. Defaults to 1.
#' @param rb The blood-to-plasma concentration ratio. Defaults to 1.
#' @param fa The fraction absorbed. Defaults to 1.
#' @param fg The fraction escaping gut metabolism. Defaults to 1.
#' @param ka The absorption rate constant in /ml. Defaults to 0.1 /ml.
#' @param oral Oral administration as Boolean. Defaults to `TRUE`.
#' @param solubility The aqueus solubility in mg/l. Defaults to Inf,
#'
#' @seealso [new_perpetrator()]
#' @seealso [read_perpetrators()]
#' @return A perpetrator object.
#' @export
#'
#' @examples
#' make_perpetrator("test", 100, 1000, 500)
make_perpetrator <- function(name, dose, imaxss, mw, type = "parent",
                             oral = TRUE, fu = 1, fumic = 1, rb = 1, fa = 1,
                             fg = 1, ka = 0.1, solubility = Inf) {
  temp <- data.frame(param = c("name", "type", "oral", "mw", "dose", "imaxss",
                               "fu", "fumic", "rb", "fa", "fg", "ka", "solubility"),
                     value = c(name, type, oral, mw, dose, imaxss, fu, fumic,
                               rb, fa, fg, ka, solubility),
                     source = rep("", 13))
  new_perpetrator(temp)
}


#' Implementation of the generic print function for perpetrator objects
#'
#' @param x The perpetrator object.
#' @param ... Further parameters.
#'
#' @export
#' @import dplyr
#' @noRd
#' @examples
#' print(examplinib_parent)
print.perpetrator <- function(x, ...) {
  cat("== DDI perpetrator object ==\n")
  x %>%
    df_to_string(colnames=F) %>%
    cat()
}



#' Generic name method
#'
#' @param obj The perpetrator object.
#'
#' @return The name of the perpetrator as character.
#' @export
#' @noRd
name <- function(obj) {
  UseMethod("name")
}


#' Name of a perpetrator compound
#'
#' @param obj The perpetrator object.
#'
#' @return The name of the perpetrator as character.
#' @export
#' @examples
#' name(examplinib_parent)
name.perpetrator <- function(obj) {
  return(obj["name", "value"])
}


#' Test if Igut of a perpetrator is limited by its solubility
#'
#' This function tests whether the solubility of a perpetrator compound is
#' limiting its intestinal concentration.
#'
#' If the 'solubility' field is `Inf` (default), or the compound solubility is
#' larger than \eqn{I_{gut}}, the function returns `FALSE`. If the solubility
#' is lower than the theoretical \eqn{I_{gut}}, i.e., lower than the dose
#' dissolved in 250 ml, the function returns `TRUE`. Note that the solubility
#' is expected in mg/l.
#'
#' @param obj The perpetrator object.
#' @return A Boolean value.
#' @seealso [key_concentrations()]
#' @export
#' @examples
#' is_igut_solubility_limited(examplinib_parent)
is_igut_solubility_limited <- function(obj) {
  oral <- as.logical(obj[which(obj$param=="oral"), "value"])
  dose <- as.num(obj[which(obj$param=="dose"), "value"])

  # total gut concentration in ng/ml
  if(oral==FALSE) {
    igut <- 0
  } else {
    igut <- dose / 250 * 1e+6
  }

  if ("solubility" %in% obj$param) {
    sol <- as.num(obj["solubility", "value"]) * 1000
    if(!is.na(sol) & igut > sol){
      return(TRUE)
    } else
      return(FALSE)
  } else {
    return(FALSE)
  }
}


#' Key perpetrator concentrations
#'
#' This function calculates the relevant perpetrator concentrations in
#' \eqn{\mu M} (default) or ng/ml for a DDI perpetrator compound.
#'
#' @details
#' ## Gut concentration
#'
#' \deqn{I_{gut} = \frac{D} {250}}
#'
#' If the above exceeds the aqueous solubility of the drug, \eqn{I_{gut}} is set
#' to its solubility.
#'
#' ## Unbound systemic concentration
#'
#' \deqn{I_{max,ss,u}=I_{max,ss} * f_u}
#'
#' ## Unbound hepatic inlet concentration
#'
#' For orally administered (parent) compounds, the hepatic inlet concentration
#' is the systemic concentration plus a portal term:
#'
#' \deqn{portal\ term = D*\frac{F_a*F_g*k_a}{Q_h*R_B}*1000\ ng/ml}
#'
#' where \eqn{D} is the administered dose in mg, \eqn{F_a} the fraction absorbed
#' after oral administration, \eqn{F_g} the fraction available after gut
#' metabolism, \eqn{k_a} the absorption rate, \eqn{Q_h} the hepatic blood flow
#' and \eqn{R_B} the blood-to-plasma ratio.
#'
#' The relevant hepatic inlet concentration (\eqn{I_{max,inlet,u}}, also called
#' \eqn{I_h} in the mechanistic static modeling equations) concentration is the
#' sum of the maximal systemic plasma concentration and the portal contribution:
#'
#' \deqn{I_{max,inlet,u}=(I_{max,ss} + portal\ term) * f_u}
#'
#' ## Unbound enteric concentration
#'
#' For orally administered (parent) compounds, the villous concentration in the
#' gut (\eqn{I_{enteric}}, also called \eqn{I_g} in the mechanistic static
#' modeling equations) is calculated as:
#'
#' \deqn{I_{enteric,u} = D * \frac{F_a*k_a}{Q_{ent}} *1000\ ng/ml}
#'
#' where \eqn{F_a} is the fraction absorbed after oral administration, \eqn{k_a}
#' the absorption rate, \eqn{Q_{ent}} the enteric villous blood flow and
#' \eqn{R_B} the blood-to-plasma distribution ratio of the compound.
#'
#' Note that as per the FDA guideline (refer to FDA, 2020, Fig. 7, and
#' Rostami-Hodjegan and Tucker, 2004) the blood-to-plasma ratio and the plasma
#' binding of the drug are ignored.
#'
#' @param qh Hepatic blood flow in l/min, defaults to 1.616 l/min.
#' @param qent Enteric blood flow in l/min, defaults to 0.3 l/min = 18 l/h.
#' @param molar Switch to select output in molar concentrations.
#' @param obj A perpetrator object.
#' @seealso [conc_table()]
#' @seealso [property_table()]
#' @return Key perpetrator concentrations in uM (default) or ng/ml as a named
#' vector.
#' @export
#' @examples
#' key_concentrations(examplinib_parent)
#' key_concentrations(examplinib_parent, molar = FALSE)
#' key_concentrations(examplinib_metabolite)
#'
key_concentrations <- function(obj, qh=1.616, qent=18/60, molar=TRUE) {
  # type <- obj["type", "value"]
  oral <- as.logical(obj["oral", "value"])
  dose <- as.num(obj["dose", "value"])
  mw <- as.num(obj["mw", "value"])
  fa <- as.num(obj["fa", "value"])
  fu <- as.num(obj["fu", "value"])
  fg <- as.num(obj["fg", "value"])
  ka <- as.num(obj["ka", "value"])
  rb <- as.num(obj["rb", "value"])
  imaxss <- as.num(obj["imaxss", "value"])

  # total gut concentration in ng/ml
  if(oral == FALSE) {
    igut <- 0
  } else {
    igut <- dose / 250 * 1e+6
  }

  if(is_igut_solubility_limited(obj)) {
    igut <- as.num(obj["solubility", "value"]) * 1000
    message(paste0("Caution: Igut for ", obj["name", "value"] ,
                   " is limited to its solubility of ",
                   obj["solubility", "value"], " mg/l!"))
  }

  # total portal contribution to hepatic inlet concentration
  if(oral == FALSE) {
    portal_term <- 0
  } else {
    portal_term <- dose * fa * fg * ka / qh / rb * 1000
  }

  # unbound systemic concentration in ng/ml
  imaxssu <- imaxss * fu

  # unbound hepatic inlet concentration in ng/ml
  imaxinletu <- (imaxss + portal_term) * fu

  # unbound intestinal concentration in ng/ml
  if(oral == FALSE) {
    imaxintestu <- imaxssu
  } else {
    imaxintestu <- dose * fa * ka / qent / rb * fu * 1000
  }

  # output vector
  temp <- c(igut=igut, imaxssu=imaxssu, imaxinletu=imaxinletu,
            imaxintestu=imaxintestu)

  if(molar) {
    temp <- temp/mw
  }
  return(temp)
}


#' Table of key perpetrator concentrations
#'
#' This function generates a markdown-formatted table of the key perpetrator
#' concentrations used for the assessment of the DDI perpetrator potential. See
#' [key_concentrations()] for details on their calculation.
#'
#' @param perp The perpetrator object or a list of perpetrator objects.
#'
#' @return A markdown-formatted table.
#' @export
#' @seealso [key_concentrations()]
#' @examples
#' conc_table(examplinib_parent)
#' conc_table(examplinib_compounds)
conc_table <- function(perp) {
  UseMethod("conc_table")
}


#' Key perpetrator concentrations
#'
#' This function generates a markdown-formatted table of the key perpetrator
#' concentrations used for the assessment of the DDI perpetrator potential. See
#' [key_concentrations()] for details on their calculation.
#'
#' @param perp The perpetrator object.
#' @return A markdown-formatted table.
#' @export
#' @seealso [key_concentrations()]
#' @examples
#' conc_table(examplinib_parent)
conc_table.perpetrator <- function(perp) {
  sol_limit <- is_igut_solubility_limited(perp)

  name <- perp["name", "value"]
  temp <- data.frame(
    parameter = c("$I_{gut}$", "$I_{max,ss,u}$", "$I_{max,inlet,u}$",
                  "$I_{max,intestinal,u}$"),
    mass_conc = format(key_concentrations(perp, molar=F), scientific=F,
                       digits=3),
    molar_conc = format(key_concentrations(perp, molar=T), scientific=F,
                        digits=3)
  )

  if (sol_limit) {
    temp["igut", "parameter"] <- "$I_{gut}$ *"
  }

  colnames(temp) <- c("parameter", "value (ng/ml)", "value (uM)")
  rownames(temp) <- NULL

  out <- knitr::kable(
    temp,
    caption = paste("Key perpetrator concentrations for", name))
  return(out)
}


#' Key perpetrator concentrations
#'
#' This function generates a list of markdown-formatted tables of the key
#' concentrations used for the assessment of the DDI perpetrator potential. See
#' [key_concentrations()] for details on the calculation of the concentrations.
#'
#' @param perp A list of perpetrator objects.
#'
#' @return A markdown-formatted table.
#' @export
#' @seealso [key_concentrations()]
#' @examples
#' conc_table(examplinib_compounds)
conc_table.list <- function(perp) {
  lapply(perp, conc_table)
}


#' Generic function to display compound properties
#'
#' @param obj The object (compound object or list thereof)
#' @export
#' @seealso [property_table.perpetrator()]
#' @seealso [property_table.list()]
property_table <- function(obj) {
  UseMethod("property_table")
}


#' Table of perpetrator drug properties
#'
#' This function generates a markdown-formatted table of the key properties of
#' a perpetrator object.
#'
#' @param obj The perpetrator object.
#'
#' @return A markdown-formatted table.
#' @export
#' @examples
#' property_table(examplinib_parent)
#' property_table(examplinib_metabolite)
#'
property_table.perpetrator <- function(obj){
  labels <- data.frame(
    param=c("mw", "dose", "solubility", "imaxss", "fu", "fumic", "rb", "fa",
            "fg", "ka"),
    parameter =c("$MW$ (g/mol)", "$dose$ (mg)", "$solubility$ (mg/l)",
                 "$C_{max,ss}$ (ng/ml)",
                 "$f_u$", "$f_{u,mic}$", "$R_B$", "$F_a$", "$F_g$",
                 "$k_a$ (1/min)")
  )

  if(as.logical(obj[which(obj$param=="oral"), "value"]) == FALSE) {
    obj <- obj %>%
      filter(!(param %in% c("fa", "fg", "ka")))
  }

  out <- obj %>%
    as.data.frame() %>%
    dplyr::slice(3:n()) %>%
    dplyr::left_join(labels, by="param") %>%
    dplyr::select(parameter, value, source) %>%
    knitr::kable(caption=paste("Compound parameters for", obj["name", "value"]))

  return(out)
}


#' Perpetrator drug properties for a list of compounds
#'
#' This function generates a list of markdown-formatted tables of the key
#' properties of a list of perpetrator objects.
#'
#' @param obj A list of perpetrator objects.
#' @return A list of markdown-formatted tables.
#' @export
#' @examples
#' property_table(examplinib_compounds)
property_table.list <- function(obj) {
  for(i in obj) {
    print(property_table(i))
  }
}







