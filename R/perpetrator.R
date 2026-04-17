#' Compound class definition
#'
#' @slot name character.
#' @slot oral logical.
#' @slot mw numeric.
#' @slot dose numeric.
#' @slot imaxss numeric.
#' @slot fu numeric.
#' @slot fumic numeric.
#' @slot rb numeric.
#' @slot fa numeric.
#' @slot fg numeric.
#' @slot ka numeric.
#' @slot solubility numeric.
#' @slot source character
#'
#' @returns A compound class object.
#' @export
setClass(
  Class = "compound",
  representation(
    name = "character",
    oral = "logical",
    mw = "numeric",
    dose = "numeric",
    imaxss = "numeric",
    fu = "numeric",
    fumic = "numeric",
    rb = "numeric",
    fa = "numeric",
    fg = "numeric",
    ka = "numeric",
    solubility = "numeric",
    source = "character"
  ),
  prototype(
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
  )
)


#' Validity check for compound class
#'
#' @param object
#'
#' @returns
#' @noRd
setValidity("compound", function(object) {
  validate_argument(object@name, "character", allow_empty = TRUE)
  validate_argument(object@oral, "logical")
  validate_argument(object@mw, "numeric", expect_positive = TRUE)
  validate_argument(object@dose, "numeric", expect_positive = TRUE)
  validate_argument(object@imaxss, "numeric", expect_positive = TRUE)
  validate_fraction(object@fu)
  validate_fraction(object@fumic)
  validate_fraction(object@rb)
  validate_fraction(object@fa)
  validate_fraction(object@fg)
  validate_argument(object@ka, "numeric", expect_positive = TRUE)
  validate_argument(object@solubility, "numeric", expect_positive = TRUE)
  validate_argument(object@source, "character", allow_multiple = TRUE, allow_empty = TRUE)

  unexpected_source <- setdiff(names(object@source), slotNames(object))
  if (length(unexpected_source) > 0)
    return(paste0("unxpected source: ", nice_enumeration(unexpected_source)))

  TRUE
})


#' Compound class constructor
#'
#' @param name Character.
#' @param oral Oral administration, as logical.
#' @param mw Molecular weight in g/mol.
#' @param dose Administered dose in mg.
#' @param imaxss Cmax at steady state in ng/ml.
#' @param fu Fraction unbound.
#' @param fumic Fraction unbound in microsomes.
#' @param rb blood cell distribution.
#' @param fa Fraction absorbed.
#' @param fg Fraction escaping gut metabolism.
#' @param ka Absorption constant in 1/min.
#' @param solubility Aqueous solubility in mg/l.
#' @param source Source information for parameters as named character vector
#'
#' @returns
#' @expor
#' @examples
#' compound(
#' "examplinib", TRUE, 492.6, 450, 3530, fu = 0.023, fa = 0.81, ka = .00267,
#' source = c(dose = "clinical dose", imaxss = "study 001", fu = "study 002")
#' )
#'
compound <- function(
    name,
    oral,
    mw,
    dose,
    imaxss,
    fu = 1,
    fumic = 1,
    rb = 1,
    fa = 1,
    fg = 1,
    ka = 0.1,
    solubility = Inf,
    source = character(0)
) {
  new(
    "compound",
    name = name,
    oral = oral,
    mw = mw,
    dose = dose,
    imaxss = imaxss,
    fu = fu,
    fumic = fumic,
    rb = rb,
    fa = fa,
    fg = fg,
    ka = ka,
    solubility = solubility,
    source = source)
}


#' Print method for compound objects.
#'
#' @param compound Compound object.
#'
#' @returns Nothing.
#' @export
setMethod(
  "show", "compound",
  function(object) {
    line <- paste0(rep("-", 5), collapse="")
    cat(paste0(line, " DDI perpetrator object ", line, "\n"))

    out <- tibble::tribble(
           ~name,   ~unit,
          "name",      "",
          "oral",      "",
            "mw", "g/mol",
          "dose",    "mg",
        "imaxss", "ng/ml",
            "fu",      "",
         "fumic",      "",
            "rb",      "",
            "fa",      "",
            "fg",      "",
            "ka",  "/min",
    "solubility", "mg/ml"
    )

    out$value <- unlist(lapply(out$name, function(x) slot(object, x)))

    out$source <- unlist(lapply(
      out$name,
      function(x) {
        ifelse(
          !is.na(slot(object, "source")[x]),
          paste0("(", slot(object, "source")[x], ")"),
          "")
      }
    ))

    out <- out |>
      mutate(value = paste0(value, " ", unit)) |>
      select(-unit) |>
      df_to_string(colnames = FALSE)

    cat(df_to_string(out, colnames = FALSE))
  }
)



#' Maximal gut concentration
#'
#' @param x A perpetrator compound object.
#' @param molar Output in molar concentration (uM). Defaults to FALSE (ng/nl).
#'
#' @returns Numeric.
#' @export
#'
igut <- function(x, molar = FALSE) {
  # total gut concentration in ng/ml
  if(oral == FALSE) {
    igut <- 0
  } else {
    igut <- x@dose / 250 * 1e+6
  }

  sol <- x@solubility * 1000
  if (!is.na(sol) &igut > sol) {
    igut <- sol
    warning(paste0("Caution: Igut is limited by solubility!"))
  }

  ifelse(molar, igut / x@mw, igut)
}


#' Cmax,ss,u
#'
#' @param x A perpetrator compound object.
#' @param molar Output in molar concentration (uM). Defaults to FALSE (ng/nl).
#'
#' @returns Numeric.
#' @export
#'
imaxssu <- function(x, molar = FALSE) {
  validate_perpetrator(x)
  out <- x@imaxss * x@fu
  ifelse(molar, out / x@mw, out)
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
  if(x@oral == FALSE) {
    0
  } else {
    x@dose * x@fa * x@fg * x@ka / qh / x@rb * 1000
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
imaxinletu <- function(x, qh = 1.616, molar = FALSE) {
  validate_perpetrator(x)
  out <- x@imaxss + portal_term(x, qh) * x@fu
  ifelse(molar, out / x@mw, out)
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
imaxintest <- function(x, qent = 18/60, molar = FALSE) {
  validate_perpetrator(x)
  if(x@oral == FALSE) {
    out <- imaxssu(x)
  } else {
    out <- x@dose * x@fa * x@ka / qent * 1000
  }

  ifelse(molar, out / x@mw, out)
}


setGeneric("kc", function(x, qh = 1.616, qent = 18/60, molar = TRUE) standardGeneric("kc"))


#' Key perpetrator concentrations
#'
#' This function returns the relevant perpetrator concentrations in
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
#' @param x A perpetrator compound object.
#' @param qh Hepatic blood flow in l/min, defaults to 1.616 l/min.
#' @param qent Enteric blood flow in l/min, defaults to 0.3 l/min = 18 l/h.
#' @param molar Output in molar concentration (uM, default).
#'
#' @returns
#' @export
#'
#' @examples
setMethod("kc", "compound", function(x, qh = 1.616, qent = 18/60, molar = TRUE) {
  validate_perpetrator(x)
  c(
    igut = igut(x, molar = molar),
    imaxssu = imaxssu(x, molar = molar),
    imaxinletu = imaxinletu(x, qh = qh, molar = molar),
    imaxintest = imaxintest(x, qent = qent, molar = molar)
  )
})



setGeneric("prop", function(x) standardGeneric("prop"))


#' Title
#'
#' @param compound
#'
#' @returns
#' @export
#'
#' @examples
setMethod("prop", "compound", function(x) {
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
    mutate(value = sapply(
      c(x@oral, x@mw, x@dose, x@solubility, x@imaxss, x@fu, x@fumic, x@rb,
        x@fa, x@fg, x@ka),
      as.character)) |>
    select(-unit)

    out$source <- unlist(lapply(
      out$param,
      function(i) {
        ifelse(
          !is.na(slot(x, "source")[i]),
          slot(x, "source")[i],
          ""
        )
      }
    ))

  if (x@oral == FALSE) {
    out <- out %>%
      filter(!(param %in% c("fa", "fg", "ka")))
  }

  out <- out %>%
    dplyr::select(parameter, value, source) %>%
    knitr::kable(caption = paste("Compound parameters for", x@name))

  out
})




setGeneric("conctbl", function(x, round = 2) standardGeneric("conctbl"))

setMethod("conctbl", "compound", function(x, round = 2) {
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
  caption <- paste0("Key perpetrator concentrations for ", x@name)

  knitr::kable(temp, col.names = col_names, caption = caption)
})


