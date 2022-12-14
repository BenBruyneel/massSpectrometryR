% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chemistry.R
\name{aminoAcidResidues}
\alias{aminoAcidResidues}
\title{Generates a pre-defined object which contains info on 'normal' amino
 acid residues}
\usage{
aminoAcidResidues()
}
\value{
a R6 object of class 'chemicals'
}
\description{
Generates a pre-defined object which contains info on 'normal' amino
 acid residues
}
\note{
The formulas in the object are amino acid residues as they are present in proteins.
 To get the actual formula of the amino acid in its 'free' form, add c(H=2, O=1) (water)

this object is used in all protein calculations in this package
}
\examples{
print(aminoAcidResidues())

}
