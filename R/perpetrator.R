#' Perpetrator class definition
#'
#' @param name character.
#' @param oral logical.
#' @param mw numeric.
#' @param dose numeric.
#' @param imaxss numeric.
#' @param fu numeric.
#' @param fumic numeric.
#' @param rb numeric.
#' @param fa numeric.
#' @param fg numeric.
#' @param ka numeric.
#' @param solubility numeric.
#' @param source character
#'
#' @returns A perpetrator class object.
#' @export
#' @examples
#' perpetrator(
#'   name = "examplinib",
#'   oral = TRUE,
#'   mw = 492.6,
#'   dose = 450,
#'   imaxss = 3530,
#'   fu = 0.023,
#'   fumic = 1,
#'   rb = 1,
#'   fa = 0.81,
#'   fg = 1,
#'   ka = 0.00267,
#'   solubility = Inf,
#'   source = c(dose = "clinical dose", imaxss = "study 001", fu = "study 002")
#' )
perpetrator <- function(
    name = "",
    oral = TRUE,
    mw = 0,
    dose = 0,
    imaxss = 0,
    fu = 1,
    fumic = 1,
    rb = 1,
    fa = 1,
    fg = 1,
    ka = 0.1,
    solubility = Inf,
    source = character(0)
){
  x <- structure(
    list(
      name = name, oral = as.logical(oral), mw = mw, dose = dose, imaxss = imaxss,
      fu = fu, fumic = fumic, rb = rb, fa = fa, fg = fg, ka = ka,
      solubility = solubility, source = source
    ),
    class = c("perpetrator", "list")
  )

  validate_perpetrator(x)

  return(x)
}


# print <- function(obj) {
#   UseMethod("print")
# }



#' Print perpetrator object as markdown
#'
#' @param x The perpetrator object.
#' @param ... Further arguments.
#'
#' @returns A markdown-formatted table.
#' @import dplyr
#' @exportS3Method ddir::print
print.perpetrator <- function(x, ...) {
  out <- tibble::tribble(
           ~param,             ~parameter,   ~unit,
           "oral",                 "oral",      "",
             "mw",         "$MW$ (g/mol)", "g/mol",
           "dose",          "$dose$ (mg)",    "mg",
     "solubility",  "$solubility$ (mg/l)",  "mg/l",
         "imaxss", "$C_{max,ss}$ (ng/ml)", "ng/ml",
             "fu",                "$f_u$",      "",
          "fumic",          "$f_{u,mic}$",      "",
             "rb",                "$R_B$",      "",
             "fa",                "$F_a$",      "",
             "fg",                "$F_g$",      "",
             "ka",        "$k_a$ (1/min)",  "/min"
     ) |>
    dplyr::mutate(value = sapply(
      c(x$oral, x$mw, x$dose, x$solubility, x$imaxss, x$fu, x$fumic, x$rb,
        x$fa, x$fg, x$ka),
      as.character)) |>
    dplyr::select(-"unit")

  out$source <- unlist(lapply(
    out$param,
    function(i) {
      ifelse(
        !is.na(x$source[i]),
        x$source[i],
        ""
      )
    }
  ))

  if (isTRUE(getOption("knitr.in.progress"))) {
    if (x$oral == FALSE) {
      out <- out |>
        dplyr::filter(!(param %in% c("fa", "fg", "ka")))
    }
    out |>
      dplyr::select(parameter, value, source) |>
      knitr::kable(caption = paste("Perpetrator compound parameters for", x$name)) |>
      print()
  } else {
    cat(paste(hline(), "DDI precipitant", x$name, hline(), "\n"))
    dplyr::select(out, -"parameter") |>
      df_to_string() |>
      cat()
  }
}


#' Maximal gut concentration
#'
#' @param x A perpetrator compound object.
#' @param molar Output in molar concentration (uM). Defaults to FALSE (ng/nl).
#'
#' @returns Numeric.
#' @export
#'
igut <- function(x, molar = FALSE) {
  validate_perpetrator(x)

  # total gut concentration in ng/ml
  oral <- x$oral
  if(oral == FALSE) {
    igut <- 0
  } else {
    igut <- x$dose / 250 * 1e+6
  }

  sol <- x$solubility * 1000
  if (!is.na(sol) &igut > sol) {
    igut <- sol
    warning(paste0("Caution: Igut is limited by solubility!"))
  }

  ifelse(molar, igut / x$mw, igut)
}


#' Cmax,ss,u
#'
#' @param x A perpetrator compound object.
#' @param molar Output in molar concentration (uM). Defaults to FALSE (ng/nl).
#'
#' @returns Numeric.
#' @export
#'
imaxssu <- function(x, molar = TRUE) {
  # input validation
  validate_perpetrator(x)

  fu = x$fu
  out <- x$imaxss * fu
  ifelse(molar, out / x$mw, out)
}



#' Portal contribution to hepatic inlet concentration
#'
#' @param x A perpetrator compound object.
#' @param qh Hepatic blood flow in l/min, defaults to 1.616 l/min.
#'
#' @returns Portal term in ng/ml
#' @noRd
#'
portal_term <- function(x, qh = 1.616) {
  validate_perpetrator(x)
  if(x$oral == FALSE) {
    0
  } else {
    x$dose * x$fa * x$fg * x$ka / qh / x$rb * 1000
  }
}


#' Cmax,u at the hepatic inlet
#'
#' @param x A perpetrator compound object.
#' @param qh Hepatic blood flow in l/min, defaults to 1.616 l/min.
#' @param molar Output in molar concentration (uM). Defaults to FALSE (ng/nl).
#'
#' @returns Numeric.
#' @export
#'
imaxinletu <- function(x, qh = 1.616, molar = TRUE) {
  validate_perpetrator(x)

  out <- (x$imaxss + portal_term(x, qh)) * x$fu
  ifelse(molar, out / x$mw, out)
}


#' Intestinal Cmax,u
#'
#' @param x A perpetrator compound object.
#' @param qent Enteric blood flow in l/min, defaults to 0.3 l/min = 18 l/h.
#' @param molar Output in molar concentration (uM). Defaults to FALSE (ng/nl).
#'
#' @returns Numeric.
#' @export
#'
imaxintest <- function(x, qent = 18/60, molar = TRUE) {
  validate_perpetrator(x)

  if(x$oral == FALSE) {
    out <- imaxssu(x, molar = molar)
  } else {
    out <- x$dose * x$fa * x$ka / qent * 1000
    ifelse(molar, out / x$mw, out)
  }

}


#' Key perpetrator concentrations
#'
#' Print a markdown-formatted table of the relevant perpetrator
#' concentrations in \eqn{\mu M} and ng/ml for a DDI perpetrator compound.
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
#' ## Enteric (villous) concentration
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
#' binding of the drug are not applicable.
#'
#' @param x DDI perpetrator object.
#' @param round Number of decimal places.
#' @param qh Hepatic blood flow in L/min.
#' @param qent Enteric blood flow in L/min.
#'
#' @returns Markdown-formatted text.
#' @export
#' @examples
#' key_conc_table(examplinib)
#'
key_conc_table <- function(
    x,
    round = 2,
    qh = 1.616,
    qent = 18/60) {
  # input validation
  validate_perpetrator(x)

  # calculations
  kc <- function(x, molar = TRUE) {
    c(
      igut = igut(x, molar = molar),
      imaxssu = imaxssu(x, molar = molar),
      imaxinletu = imaxinletu(x, qh = qh, molar = molar),
      imaxintest = imaxintest(x, qent = qent, molar = molar)
    )
  }

  temp <- tibble::tribble(
    ~parameter,
    "$I_{gut}$",
    "$I_{max,ss,u}$",
    "$I_{max,inlet,u}$",
    "$I_{max,intestinal}$"
  ) |>
    mutate(mass_conc = round(kc(x, molar = FALSE), round)) |>
    mutate(molar_conc = round(kc(x, molar = TRUE), round))

  col_names <- c("parameter", "value ($ng/ml$)", "value ($\\mu M$)")
  rownames(temp) <- NULL
  caption <- paste0("Key perpetrator concentrations for ", x$name)

  knitr::kable(temp, col.names = col_names, caption = caption)
}


