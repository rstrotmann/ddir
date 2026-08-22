#' Horizontal line using unicode characters
#'
#' @param length Line length.
#'
#' @returns A character string.
#' @noRd
hline <- function(length = 8) {
  paste0(rep("\U2500", length), collapse = "")
}


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
#' @param label_empty annotate empty data set, as logical.
#'
#' @return The data frame representation as character.
#' @import stringr
#' @import utils
#' @noRd
df_to_string <- function(df, indent = "", n = NULL, colnames = TRUE,
                         label_empty = TRUE){
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

  if (nrow(df) > 0) {
    for(i in 1:nrow(df)){
      out <- paste(out, render.line(df[i,]), sep="\n")
    }
  } else {
    if (label_empty == TRUE)
      out <- paste(out, "(empty data set)", sep="\n")
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

#' Nice enumeration of multiple strings
#'
#' @param items Items to enumerate as character.
#' @param conjunction The conjunction between the last and penultmate items.
#'
#' @return Enumeration as character.
#' @noRd
#'
#' @examples
#' nice_enumeration("A")
#' nice_enumeration(c("A", "B"))
#' nice_enumeration(c("A", "B", "C"))
#' nice_enumeration(c("A", "B", "C"), conjunction = "or")
nice_enumeration <- function(items, conjunction = "and") {
  if (length(items) == 1) {
    return(items[[1]])
  }
  if (length(items) > 1) {
    # items <- Filter(function(x) x != "", items)
    return(paste(
      paste(items[1:length(items) - 1], collapse = ", "), conjunction,
      items[length(items)]
    ))
  }
}


#' Make column names in markdown format
#'
#' @param obj A data frame.
#'
#' @returns A character vector.
#' @noRd
make_labels <- function(obj) {
  translation <- tibble::tribble(
             ~name,                 ~label,
              "ki",   "$K_{i}$ ($\\mu M$)",
            "ic50",            "$IC_{50}$",
               "i",                  "$i$",
             "kiu", "$K_{i,u}$ ($\\mu M$)",
               "r",                  "$R$",
           "r_gut",            "$R_{gut}$",
        "risk_hep",       "risk (hepatic)",
     "risk_intest",    "risk (intestinal)",
          "kinact",    "$k_{inact}$ (1/h)",
            "kdeg",      "$k_{deg}$ (1/h)",
              "fu",                "$f_u$",
            "emax",            "$E_{max}$",
            "ec50",            "$EC_{50}$",
           "max_c",   "$max c$ ($\\mu M$)",
            "aucr",            "$R_{AUC}$",
           "fmcyp",           "$fm_{CYP}$",
            "fgut",            "$f_{gut}$"
     )


  as.character(
    sapply(
      names(obj),
      function(x) {
        temp <- translation[translation$name == x,]
        ifelse(nrow(temp) == 0, x, temp$label)
      }
    )
  )
}





# c("CYP", "$E_{max}$", "$max c$ ($\\mu M$)", "source",
#   "$max c/C_{max,ss,u}$", "risk", "notes")
