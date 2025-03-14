% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary_statements.R
\name{ugt_inhibition_risk_summary}
\alias{ugt_inhibition_risk_summary}
\title{Summary statement for UGT inhibition}
\usage{
ugt_inhibition_risk_summary(perp, ugt_inh)
}
\arguments{
\item{perp}{The perpetrator object.}

\item{ugt_inh}{UGT inhibition data as data frame, The following fields are
expected:
\itemize{
\item 'name' The name of the perpetrator compound.
\item 'ugt' The UGT enzyme as (upper case) character.
\item 'ic50' The \eqn{IC_{50}} in \eqn{\mu M} as numeric.
\item 'source' Optional source information as character.
}}
}
\value{
Markdown-formatted text output.
}
\description{
Summary statement for UGT inhibition
}
\examples{
ugt_inhibition_risk_summary(
  examplinib_parent, examplinib_ugt_inhibition_data)
ugt_inhibition_risk_summary(
  examplinib_compounds, examplinib_ugt_inhibition_data)
}
