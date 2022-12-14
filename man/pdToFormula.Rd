% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chemistry.R
\name{pdToFormula}
\alias{pdToFormula}
\title{translates a proteome Discoverer (Thermo Scientific) elements formula string
 to a formula as used by this package}
\usage{
pdToFormula(pdFormula)
}
\arguments{
\item{pdFormula}{a character vector. Formula in a format as used by proteome
discoverer software}
}
\value{
formula of format c(H=2, O=1)
}
\description{
translates a proteome Discoverer (Thermo Scientific) elements formula string
 to a formula as used by this package
}
\examples{
glucose <- pdToFormula("C(6) H(12) O(6)")
glucose
water <- pdToFormula("H(2) O")
water
}
