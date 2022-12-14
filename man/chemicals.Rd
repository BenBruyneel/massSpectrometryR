% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chemistry.R
\name{chemicals}
\alias{chemicals}
\title{R6 Class representing a set of chemicals}
\description{
Every chemical inside the object has a name, letter, short and a formula.
 The first 3 can be any length of string (though the letter and short field
  should be maximum length (nchar) 1 and 2-4 respectively). Formula should be
  in the form of a named numeric with the names representing elements and the
  values themselves being the number of atoms of that element,
  eg c(C = 3, H =  5, N = 1, O = 1, S = 0)
}
\details{
Note: this class is meant to be used for classes of compounds, eg amino acids

Also: his class is  meant as a base class to be expanded via inheritance

Warning: all chemicals inside this object should be unique (names, letters &
 shorts)
}
\examples{
estrogens <- chemicals$new(letters = c("1","2","3","4"),
                           shorts = c("E1","E2","E3","E4"),
                           names = c("Estrone","Estradiol",
                                     "Estriol","Estetrol"),
                           formulas = list(c(C=18, H=22, O=2),
                                           c(C=18, H=24, O=2),
                                           c(C=18, H=24, O=3),
                                           c(C=18, H=24, O=4)))

}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{number}}{retrieve the number of compounds present in the object, read only}

\item{\code{letters}}{to access the letters of the compounds in the object}

\item{\code{names}}{to access the names of the compounds in the object}

\item{\code{shorts}}{to access the shorts of the compounds in the object}

\item{\code{formulas}}{to access the formulas of the compounds in the object}

\item{\code{table}}{retrieves all info on the compounds in data.frame format, read only}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-chemicals-new}{\code{chemicals$new()}}
\item \href{#method-chemicals-print}{\code{chemicals$print()}}
\item \href{#method-chemicals-getFormula}{\code{chemicals$getFormula()}}
\item \href{#method-chemicals-clone}{\code{chemicals$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-chemicals-new"></a>}}
\if{latex}{\out{\hypertarget{method-chemicals-new}{}}}
\subsection{Method \code{new()}}{
Create a new chemicals object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{chemicals$new(letters, shorts, names, formulas)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{letters}}{character vector specifying the letters (or numbers or whatever) for the chemicals.
In case of amino acids it should be eg "A" for Alanine, "G" for Glycine, etc etc}

\item{\code{shorts}}{character vector specifying the short names for the chemicals, eg Ala for Alanine}

\item{\code{names}}{character vector specifying the names of the chemicals}

\item{\code{formulas}}{list of named numeric vectors specifying the formulas of the chemicals,
eg c(C = 6, H = 12, N = 4, O = 1, S = 0) for Arginine}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a new 'chemical' object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-chemicals-print"></a>}}
\if{latex}{\out{\hypertarget{method-chemicals-print}{}}}
\subsection{Method \code{print()}}{
For printing purposes: prints a table of the chemicals with columns letter, name & short
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{chemicals$print(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{no arguments, the function takes care of printing}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-chemicals-getFormula"></a>}}
\if{latex}{\out{\hypertarget{method-chemicals-getFormula}{}}}
\subsection{Method \code{getFormula()}}{
Retrieves the formula of one of the compounds in the object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{chemicals$getFormula(which1)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{which1}}{specifies which chemical should be retrieved. Which the number (row number in the chemicals
table), or the name, letter or short as a character vector. The way this is set up, it doesn't matter
whether capital or non-capital letters are used, since all is converted to upper case before comparing
with what's in the chemical table}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a formula in the shape of a named numeric vector, eg c(C = 6, H = 12, N = 4, O = 1, S = 0)
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-chemicals-clone"></a>}}
\if{latex}{\out{\hypertarget{method-chemicals-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{chemicals$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
