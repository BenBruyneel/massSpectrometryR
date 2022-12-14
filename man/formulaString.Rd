% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chemistry.R
\name{formulaString}
\alias{formulaString}
\title{Translates regular formula format into a character vector, eg
 C6H12O6}
\usage{
formulaString(formula)
}
\arguments{
\item{formula}{named numeric vector, example c(O = 2, C = 1)}
}
\value{
character vector
}
\description{
Translates regular formula format into a character vector, eg
 C6H12O6
}
\examples{
formulaString(c(C=6, H=12, O=6))
formulaString(c(H=3,O=4,P=1))
}
