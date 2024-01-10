#' Render data frame to character.
#'
#' This function renders a data frame into a string object.
#'
#' @param df The datavframe.
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
#' @param perps perpetrator objects as a list.
#' @import lifecycle
#' @return The output as string.
#' @export
names_string <- function(perps) {
  lifecycle::deprecate_warn("0.8.1", "names_string()", "compound_names_string()")
  return(paste(lapply(perps, function(p) {p[(p$param=="name"), "value"]}), collapse=", "))
}


#' Compound names as single string
#'
#' This function returns a single string with the comma-separated names of the
#' compounds included in the list of compunds.
#'
#' @param compounds A list of perpetrator objects
#'
#' @return A character string.
#' @export
#'
#' @examples
#' compound_names_string(examplinib_compounds)
compound_names_string <- function(compounds) {
  return(paste(lapply(compounds, function(p) {
    p[(p$param=="name"), "value"]}),
    collapse=", "))
}


#' Perpetrator object constructor
#'
#' @details
#'
#' The input is a data frame with the columns 'name', 'param', 'value' and
#' 'source'.
#' In all rows, 'name' must be the name of the compound, 'param' defines the
#' parameter, 'value' is the respective value and 'source' is (optional) source
#' information, e.g., the code of the study from which the value is taken.
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
#' | ka | optional | absorption rate constant, default is 0.1 /ml |
#' | solubility | optional | solubility of the compound in mg/l, defaults to `Inf` |
#'
#' The following example is a valid input data frame:
#'
#' |       name |      param |      value |        source |
#' | --- | --- | --- | --- |
#' | examplinib |       name | examplinib |               |
#' | examplinib |       type |     parent |               |
#' | examplinib |         mw |      492.6 |               |
#' | examplinib |       dose |        450 | clinical dose |
#' | examplinib |     imaxss |       3530 |     study 001 |
#' | examplinib |         fu |      0.023 |     study 002 |
#' | examplinib |      fumic |          1 |       default |
#' | examplinib |         rb |          1 |     study 003 |
#' | examplinib |         fa |       0.81 |     study 003 |
#' | examplinib |         fg |          1 |       default |
#' | examplinib |         ka |    0.00267 |       unknown |
#' | examplinib | solubility |        Inf |       default |
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
    "type",   "parent",
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

  df %>%
    assertr::verify(assertr::has_all_names("param", "value", "source")) %>%
    assertr::verify(c("name", "mw", "dose", "imaxss") %in% .$param)

  compound_name <- df[which(df$param=="name"), "value"]

  out <- df %>%
    as.data.frame() %>%
    assertr::verify(name == compound_name) %>%
    right_join(default_values, by="param") %>%
    tidyr::fill(name) %>%
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


#' Implementation of the generic print function for perpetrator objects
#'
#' @param x The perpetrator object.
#' @param ... Further parameters.
#'
#' @export
#' @import dplyr
#' @examples
#' print(examplinib_parent)
#'
print.perpetrator <- function(x, ...) {
  cat("== DDI perpetrator object ==\n")
  x %>%
    dplyr::select(-name) %>%
    df_to_string(colnames=F) %>%
    cat()
}



#' Generic name method
#'
#' @param obj The perpetrator object.
#'
#' @return The name of the perpetrator as character.
#' @export
name <- function(obj) {
  UseMethod("name")
}


#' Name of a perpetrator
#'
#' @param obj The perpetrator object.
#'
#' @return The name of the perpetrator as character.
#' @export
#' @examples
#' name(examplinib_parent)
#'
name.perpetrator <- function(obj) {
  return(obj["name", "value"])
}


#' Test if Igut of a perpetrator is limited by its solubility
#'
#' If the solubility field is `Inf` (default), or the compound solubility is
#' larger than Igut, the function returns `FALSE`. If the solubility is lower
#' than the theoretical Igut, i.e., lower than the dose dissolved in 250 ml,
#' the function returns `TRUE`. Note that the solubility is expected in mg/l.
#'
#' @param obj The perpetrator object.
#'
#' @return A boolean value.
#' @export
#' @examples
#' is_igut_solubility_limited(examplinib_parent)
#'
is_igut_solubility_limited <- function(obj) {
  type <- obj[which(obj$param=="type"), "value"]
  dose <- as.num(obj[which(obj$param=="dose"), "value"])

  # total gut concentration in ng/ml
  if(type=="metabolite") {
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
#' @details
#'
#' ## Gut concentration
#'
#' \deqn{I_{gut} = \frac{D} {250}}
#'
#' ## Systemic concentration
#'
#' \deqn{I_{max,ss,u}=I_{max,ss} * f_u}
#'
#' ## Hepatic inlet concentration
#'
#' For orally administered (parent) compounds, the hepatic inlet concentration
#' is the sytemic concentration plus a portal term:
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
#' ## Enteric concentration
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
#' @param molar Boolean value to select output in molar concentrations.
#' @param obj A perpetrator object.
#' @seealso [conc_table()]
#' @seealso [property_table()]
#' @return Key perpetrator concentrations as a named vector.
#' @export
#' @examples
#' key_concentrations(examplinib_parent)
#' key_concentrations(examplinib_metabolite)
#'
key_concentrations <- function(obj, qh=1.616, qent=18/60, molar=TRUE) {
  type <- obj["type", "value"]
  dose <- as.num(obj["dose", "value"])
  mw <- as.num(obj["mw", "value"])
  fa <- as.num(obj["fa", "value"])
  fu <- as.num(obj["fu", "value"])
  fg <- as.num(obj["fg", "value"])
  ka <- as.num(obj["ka", "value"])
  rb <- as.num(obj["rb", "value"])
  imaxss <- as.num(obj["imaxss", "value"])

  # total gut concentration in ng/ml
  if(type=="metabolite") {
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
  if(type=="metabolite"){
    portal_term <- 0
  } else {
    portal_term <- dose * fa * fg * ka / qh / rb * 1000
  }

  # unbound systemic concentration in ng/ml
  imaxssu <- imaxss * fu

  # unbound hepatic inlet concentration in ng/ml
  imaxinletu <- (imaxss + portal_term) * fu

  # unbound intestinal concentration in ng/ml
  if(type=="metabolite") {
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


#' Generic function to display ky perpetrator concentrations
#'
#' @param perp The object (compund object or list thereof)
#'
#' @export
#' @seealso [key_concentrations()]
conc_table <- function(perp) {
  UseMethod("conc_table")
}


#' Table of key perpetrator concentrations
#'
#' @param perp The perpetrator object.
#'
#' @return Key perpetrator concentrations as a markdown-formatted table.
#' @export
#' @seealso [key_concentrations()]
#' @examples
#' conc_table(examplinib_parent)
conc_table.perpetrator <- function(perp) {
  sol_limit <- is_igut_solubility_limited(perp)

  name <- perp["name", "value"]
  temp <- data.frame(
    parameter = c("$I_{gut}$", "$I_{max,ss,u}$", "$I_{max,inlet,u}$", "$I_{max,intestinal,u}$"),
    mass_conc = format(key_concentrations(perp, molar=F), scientific=F, digits=3),
    molar_conc = format(key_concentrations(perp, molar=T), scientific=F, digits=3)
  )

  if (sol_limit) {
    temp["igut", "parameter"] <- "$I_{gut}$ *"
  }

  colnames(temp) <- c("parameter", "value (ng/ml)", "value (uM)")
  rownames(temp) <- NULL

  out <- knitr::kable(temp, caption = paste("Key perpetrator concentrations for", name))
  return(out)
}


#' Table of key perpetrator concentrations for a list of compounds
#'
#' @param perp The compounds as list of compound objects.
#'
#' @export
#' @examples
#' conc_table(examplinib_compounds)
conc_table.list <- function(perp) {
  lapply(perp, conc_table)
}


#' Generic function to display compound properties
#'
#' @param obj The object (compund object or list thereof)
#'
#' @export
property_table <- function(obj) {
  UseMethod("property_table")
}


#' Table of perpetrator properties
#'
#' @param obj The perpetrator object.
#'
#' @return Perpetrator properties as markdown-formatted table.
#' @export
#' @examples
#' property_table(examplinib_parent)
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

  out <- obj %>%
    as.data.frame() %>%
    dplyr::slice(3:n()) %>%
    dplyr::left_join(labels, by="param") %>%
    dplyr::select(parameter, value, source) %>%
    knitr::kable(caption=paste("Compound parameters for", obj["name", "value"]))

  return(out)
}


#' Compound properties for a list of compounds
#'
#' @param obj A list of perpetrator objects.
#'
#' @export
#' @examples
#' property_table(examplinib_compounds)
property_table.list <- function(obj) {
  for(i in obj) {
    print(property_table(i))
  }
}







