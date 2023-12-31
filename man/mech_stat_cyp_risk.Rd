% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cyp.R
\name{mech_stat_cyp_risk}
\alias{mech_stat_cyp_risk}
\title{CYP perpetration risk as per mechanistic-static modeling}
\usage{
mech_stat_cyp_risk(
  perp,
  cyp_inh,
  cyp_ind,
  include_induction = T,
  substr = cyp_reference_substrates
)
}
\arguments{
\item{perp}{The perpetrator object.}

\item{cyp_inh}{CYP inhibiton data as data frame.}

\item{cyp_ind}{CYP induction data as data frame.}

\item{include_induction}{Boolean value to define whether induction should be
included in the calculation (C-terms as per the FDA guideline)}

\item{substr}{The CYP reference substrates to be used as data frame.}
}
\value{
A data frame.
}
\description{
CYP perpetration risk as per mechanistic-static modeling
}
