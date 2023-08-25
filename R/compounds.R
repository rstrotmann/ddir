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
