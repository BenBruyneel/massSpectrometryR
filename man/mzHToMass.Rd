% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chemistry.R
\name{mzHToMass}
\alias{mzHToMass}
\title{calculates the mass of the molecule in an ion (M+xH)x+}
\usage{
mzHToMass(mz, charge = 1, elementsInfo = elementsMonoisotopic())
}
\arguments{
\item{mz}{numeric vector, mass to charge ratio of he ion}

\item{charge}{charge state}

\item{elementsInfo}{elements masses to be used, needs to be of class
elements, default is elementsMonoisotopic()}
}
\value{
numeric vector
}
\description{
a wrapper around \code{\link{mzToMass}} for positively charged,
 protonated ions in ESI
}
\examples{
 massToMzH(mass = 174.1117, charge = 2) |> mzHToMass(charge = 2)
 massToMzH(mass = 174.1117, charge = 1) |> mzHToMass(charge = 1)
}
