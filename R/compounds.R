#' Render data frame object to string
#'
#' This function renders a data.frame into a string similar to its representation
#'  when printed without line numbers
#'
#' @param df The data.frame to be rendered
#' @param indent A string that defines the left indentation of the rendered
#'   output.
#' @param n The number of lines to be included, or all if NULL.
#' @return The output as string.
#' @import stringr
#' @import utils
df_to_string <- function(df, indent="", n=NULL, colnames=TRUE){
  df <- as.data.frame(df)
  max.widths <- as.numeric(lapply(rbind(df, names(df)), FUN=function(x) max(sapply(as.character(x), nchar), na.rm=TRUE)))
  line = df[1,]

  render.line <- function(line){
    out <- indent
    for(i in 1:length(line)){
      out <- paste0(out, sprintf(paste0("%-", max.widths[i]+3, "s"), as.character(line[i])))
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

#' Load perpetrator file
#'
#' @param filename The full path to the compound file
#'
#' @return A list of perpetrator objects
#' @export
#' @import dplyr
#' @import utils
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
  out <- lapply(data, perpetrator)
  return(out)
}

#' @export
perpetrator <- function(df) {
  stopifnot(c("param", "value", "source") %in% colnames(df))
  rownames(df) <- df$param
  stopifnot(c("name", "type", "mw", "dose", "imaxss", "fu", "fumic", "rb",
                  "fa", "fg", "ka") %in% rownames(df))
  class(df) <- c("perpetrator", "data.frame")
  df
}

#' @export
#' @import dplyr
print.perpetrator <- function(obj) {
  cat("== DDI perpetrator object ==\n")
  obj %>%
    dplyr::select(-name) %>%
    df_to_string(colnames=F) %>%
    cat()
}


#' Name of a perpetrator
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
name.perpetrator <- function(obj) {
  return(obj["name", "value"])
}



#' Test if Igut of a perpetrator is limited by its solubility
#'
#' @param obj The perpetrator object.
#'
#' @return A boolean value.
is_igut_solubility_limited <- function(obj) {
  type <- obj["type", "value"]
  dose <- as.num(obj["dose", "value"])

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
    }
  } else {
    return(FALSE)
  }
}


#' Key perpetrator concentrations
#'
#' @param qh Hepatic blood flow in l/min.
#' @param qent Enteric blood flow in l/min.
#' @param molar Boolean value to select output in molar concentrations.
#' @param obj A perpetrator object.
#'
#' @return Key perpetrator concentrations as a named vector.
#' @export
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


#' Table of key perpetrator concentrations
#'
#' @param perp The perpetrator object.
#'
#' @return Key perpetrator concentrations as a markdown-formatted table.
#' @export
conc_table <- function(perp) {
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

  colnames(temp) <- c("parameter", "value (ng/ml)", "value (µM)")
  rownames(temp) <- NULL

  out <- knitr::kable(temp, caption = paste("Key perpetrator concentrations for", name))
  # if (sol_limit) {
  #   out <- kableExtra::add_footnote(
  #     out,
  #     label="*$I_{gut}$ is limited by the solubility of the compound",
  #     notation="none",
  #     threeparttable = TRUE)
  # }
  return(out)
}



#' Table of perpetrator properties
#'
#' @param obj The perpetrator object.
#'
#' @return Perpetrator properties as markdown-formatted table.
#' @export
property_table <- function(obj){
  labels <- data.frame(
    param=c("mw", "dose", "solubility", "imaxss", "fu", "fumic", "rb", "fa", "fg", "ka"),
    parameter =c("$MW$ (g/mol)", "$dose$ (mg)", "$solubility$ (mg/l)", "$C_{max,ss}$ (ng/ml)",
            "$f_u$", "$f_{u,mic}$", "$R_B$", "$F_a$", "$F_g$", "$k_a$ (1/min)")
  )

  out <- obj %>%
    as.data.frame() %>%
    dplyr::slice(3:n()) %>%
    dplyr::left_join(labels, by="param") %>%
    dplyr::select(parameter, value, source) %>%
    knitr::kable(caption=paste("Compound parameters for", obj["name", "value"]))

  return(out)
}


