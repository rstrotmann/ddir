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
