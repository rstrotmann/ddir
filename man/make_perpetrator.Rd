% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compounds.R
\name{make_perpetrator}
\alias{make_perpetrator}
\title{Make a perpetrator object from compound data}
\usage{
make_perpetrator(
  name,
  dose,
  imaxss,
  mw,
  type = "parent",
  oral = TRUE,
  fu = 1,
  fumic = 1,
  rb = 1,
  fa = 1,
  fg = 1,
  ka = 0.1,
  solubility = Inf
)
}
\arguments{
\item{name}{The compound name as character.}

\item{dose}{The clinical dose in mg.}

\item{imaxss}{The (total) steady-state Cmax in ng/ml.}

\item{mw}{The molar weight in g/mol.}

\item{type}{The compound type as character. Must be either 'parent' or
metabolite'.}

\item{oral}{Oral administration as Boolean. Defaults to \code{TRUE}.}

\item{fu}{The fraction unbound as numeric. Defaults to 1.}

\item{fumic}{The fraction unbound in the microsomes. Defaults to 1.}

\item{rb}{The blood-to-plasma concentration ratio. Defaults to 1.}

\item{fa}{The fraction absorbed. Defaults to 1.}

\item{fg}{The fraction escaping gut metabolism. Defaults to 1.}

\item{ka}{The absorption rate constant in /ml. Defaults to 0.1 /ml.}

\item{solubility}{The aqueus solubility in mg/l. Defaults to Inf,}
}
\value{
A perpetrator object.
}
\description{
This function creates a perpetrator object from key compound data. Note that
the 'source' field is empty. If you want to include source information, you
can either use \code{\link[=new_perpetrator]{new_perpetrator()}} which takes as input a data that may
include source information, or you can create a perpetrator object from a
text input string using \code{\link[=read_perpetrators]{read_perpetrators()}}.
}
\examples{
make_perpetrator("test", 100, 1000, 500)
}
\seealso{
\code{\link[=new_perpetrator]{new_perpetrator()}}

\code{\link[=read_perpetrators]{read_perpetrators()}}
}
