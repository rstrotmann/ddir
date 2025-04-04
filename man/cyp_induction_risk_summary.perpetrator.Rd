% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary_statements.R
\name{cyp_induction_risk_summary.perpetrator}
\alias{cyp_induction_risk_summary.perpetrator}
\title{Summary statement for CYP induction for a precipitant compound}
\usage{
\method{cyp_induction_risk_summary}{perpetrator}(perp, cyp_ind, d = 1)
}
\arguments{
\item{perp}{The perpetrator object or a list thereof.}

\item{cyp_ind}{The CYP induction data as data frame. The following fields
are expected:
\itemize{
\item 'name' The name of the perpetrator compound as character.
\item 'cyp' The CYP enzyme as (upper case) character.
\item 'emax' The \eqn{E_{max}}, i.e., the maximum induction effect determined in
vitro as numeric.
\item 'ec50' The \eqn{EC_{50}} in \eqn{\mu M} as numeric.
\item 'maxc' The maximal concentration in \eqn{\mu M} tested in the in vitro assay as
numeric.
\item 'source' Optional source information as character.
}}

\item{d}{Scaling factor, defaults to 1.}
}
\value{
Markdown-formatted text output.
}
\description{
Summary statement for CYP induction for a precipitant compound
}
