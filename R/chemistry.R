# ---- elements ----

#' R6 Class representing a set of elements
#'
#' Note: this class is meant to be used for elements.
#'
#' Warning: all elements inside this object should be unique (names & shorts,
#'  not mass)
#'
#' @description
#' Every element inside the object has a name, letter and a mass. The first 2
#'  can be any length of string, mass should be a numeric
#'
#' @examples
#' randomElements <- elements$new(shorts = c("X1","X2","X3"),
#'                                names  = c("Secret Element 1",
#'                                           "Secret Element 2",
#'                                           "Secret Element 1"),
#'                                mass   = c(301, 312, 323))
#'
#' @export
elements <- R6::R6Class("elements",
                        private = list(
                          # names: the actual names of the elements
                          names_ = NA,
                          # shorts: the short names for the elements,
                          #  eg Hg for Mercury
                          shorts_ = NA,
                          # mass: the masses for the elements
                          mass_ = NA
                        ),
                        public = list(
                          #' @description creates a new elements object
                          #'
                          #' @param shorts character vector specifying the short names for the elements, eg Hg for Mercury
                          #' @param names character vector specifying the names of the elements
                          #' @param mass numeric vector specifying the masses of the elements
                          #'
                          #' @return a new 'elements' object
                          #' @export
                          initialize = function(shorts, names, mass){
                            stopifnot(length(shorts) == length(names) &
                                        length(shorts) == length(mass))
                            private$names_ <- names
                            private$shorts_ <- shorts
                            private$mass_ <- mass
                            invisible(self)
                          },
                          #' @description adds one or more elements to the object. Elements to be added must have unique
                          #'  names and shorts.
                          #'
                          #' @param shorts character vector specifying the short names for the elements, eg Hg for Mercury
                          #' @param names character vector specifying the names of the elements
                          #' @param mass numeric vector specifying the masses of the elements
                          #'
                          #' @return nothing
                          #' @export
                          addElement = function(shorts, names, mass){
                            stopifnot(length(shorts) == length(names) &
                                        length(shorts) == length(mass))
                            if ((names %in% private$names_) |
                                (shorts %in% private$shorts_)){
                              stop("Element short or name already is in object")
                            }
                            private$names_ <- append(private$names_, names)
                            private$shorts_ <- append(private$shorts_,shorts)
                            private$mass_ <- append(private$mass_, mass)
                            invisible(self)
                          },
                          #' @description
                          #' For printing purposes: prints a table of the chemicals with columns letter, name & short
                          #' @param ... no arguments, the function takes care of printing
                          print = function(...){
                            print(self$table)
                          },
                          #' @description Retrieves the mass of one of the elements in the object
                          #' @param which1 specifies which element should be retrieved. Which number (row number in
                          #'  the chemicals table), name or short as a character vector. The way this is set up,
                          #'  it doesn't matter whether capital or non-capital letters are used, since all are
                          #'  converted to upper case before comparing with what's in the elements table
                          #' @return a numeric value
                          getMass = function(which1){
                            stopifnot(!is.null(which1) & !is.na(which1))
                            # number
                            if (which1 %in% 1:self$number){
                              return(self$mass[which1])
                            } else {
                              which1 <- toupper(which1)
                              # letter
                              temporary <- which(toupper(self$shorts) == which1)
                              if (!(identical(temporary, integer(0)))){
                                return(self$mass[temporary])
                              } else {
                                # short
                                temporary <- which(toupper(self$names) == which1)
                                if (!(identical(temporary, integer(0)))){
                                  return(self$mass[temporary])
                                } else {
                                  return(NA)
                                }
                              }
                            }
                          }
                        ),
                        active = list(
                          #' @field number retrieve the number of elements present in the object, read only
                          number = function(value){
                            if (missing(value)){
                              return(length(private$names_))
                            } else {
                              # nothing, readonly!
                            }
                          },
                          #' @field names to access the names of the elements in the object
                          names = function(value){
                            if (missing(value)){
                              return(private$names_)
                            } else {
                              private$names_ <- value
                            }
                          },
                          #' @field shorts to access the shorts of the elements in the object
                          shorts = function(value){
                            if (missing(value)){
                              return(private$shorts_)
                            } else {
                              private$shorts_ <- value
                            }
                          },
                          #' @field mass to access the masses of the elements in the object
                          mass = function(value){
                            if (missing(value)){
                              return(private$mass_)
                            } else {
                              private$mass_ <- value
                            }
                          },
                          #' @field table retrieves all info on the elements in data.frame format, read only
                          table = function(value){
                            if (missing(value)){
                              return(data.frame(
                                name = private$names_,
                                short = private$shorts_,
                                mass = private$mass_))
                            } else {
                              # nothing, readonly
                            }
                          }
                        )
)

#' @title generates a pre-defined R6 elements object which contains info on
#'  elements, mass values are mono isotopic
#' @note this object is (by default) used in all chemical calculations in this
#'  package
#' @note electron is not really meant to be used in formulas, but is needed to
#'  calculate mass & m/z of ions
#' @return an elements object
#' @export
#'
#' @examples
#' print(elementsMonoisotopic())
elementsMonoisotopic <- function(){
  return(
    elements$new(
      shorts = c("C","H","N","O","S",
                "P","Br","Cl","F","Si",
                "e","Na","K","13C","15N",
                "18O","2H"),
     names  = c("Carbon", "Hydrogen","Nitrogen","Oxygen","Sulphur",
                "Phosphorus","Bromine","Chlorine","Fluorine","Silicon",
                "electron", "Sodium","Potassium","Carbon13","Nitrogen15",
                "Oxygen18","Deuterium"),
    mass   = c(12.00000000,  1.007825032, 14.00307401, 15.99491462, 31.97207073,
               30.97376149, 78.9183379  , 34.96885271, 18.99840320, 27.97692649,
                0.00054858, 22.989769660, 38.96370690, 13.00335484, 15.00010897,
               17.9991604, 2.014101778))
  )
}

#' @title generates a pre-defined object which contains info on elements,
#'  mass values are average masses (weighted mean mass of elements based on
#'  their natural occurence)
#' @note electron is not really meant to be used in formulas, but is needed to
#'  calculate mass & m/z of ions
#' @return an elements object
#' @export
#'
#' @examples
#' print(elementsAverage())
elementsAverage <- function(){
  return(
    elements$new(
      shorts = c("C","H","N","O","S",
                 "P","Br","Cl","F","Si",
                 "e","Na","K"),
      names  = c("Carbon","Hyrdogen","Nitrogen","Oxygen","Sulphur",
                 "Phosphorus","Bromine","Chlorine","Fluorine", "Sillicon",
                 "electron","Sodium","Potassium"),
      mass   = c(12.0107,  1.0079, 14.0067, 15.9994, 32.0648,
                 30.9738, 79.9035, 35.4529, 18.9984, 28.0855,
                  0.0005, 22.9898, 39.0983))
  )
}

# ---- calculations/Formulas ----

#' @title sortFormula
#' @description sorts the elements of a formula in alphabetical order
#'  (increasing/decreasing)
#'
#' @param formula named numeric vector, example c(O = 2, C = 1)
#' @param decrease logical flag on how to sort, default = FALSE: increasing
#'
#' @return named numeric vector (formula)
#' @export
#'
#' @examples
#' glucose = c(O=6, H=12, C=6)
#' glucose
#' sortFormula(glucose)
sortFormula <- function(formula, decrease = FALSE){
  return(formula[order(names(formula), decreasing = decrease)])
}

#' @title removeZeros
#' @description removes elements that have number zero
#'
#' @param formula named numeric vector, example c(O = 2, C = 1)
#'
#' @return named numeric vector (formula)
#' @export
#'
#' @examples
#' glucose = c(O=6, H=12, C=6)
#' glucose %f+% emptyFormula()
#' removeZeros(glucose %f+% emptyFormula())
removeZeros <-function(formula){
  return(formula[formula != 0])
}

#' @title elementsInFormula
#' @description retrieves the names of the elements in a formula
#'
#' @param formula named numeric vector, example c(O = 2, C = 1)
#' @param removeZero logical flag on how to deal with elements which are zero
#'
#' @return character vector
#' @export
#'
#' @examples
#' glucose = c(C=6, H=12, O=6, S=0)
#' elementsInFormula(glucose)
#' elementsInFormula(glucose, removeZero = TRUE)
elementsInFormula <- function(formula, removeZero = FALSE){
  if (removeZero){
    formula <- formula[formula != 0]
  }
  return(names(formula))
}

#' @title Combine elements present in 2 separate formulas
#' @description especially important when formula1 contains elements that
#'  formula2 does not contain and vice versa, Note: also sorts the result
#'
#' @param formula1 named numeric vector, example c(O = 2, C = 1)
#' @param formula2 named numeric vector, example c(H = 2, S = 1)
#' @param decrease logical flag on how to sort, default = FALSE: increasing
#'
#' @return a character vector
#' @export
#'
#' @examples
#' elementsInFormulas(c(O = 2, C = 1),
#'                    c(H = 2, S = 1))
elementsInFormulas <- function(formula1, formula2, decrease = FALSE){
  e1 <- elementsInFormula(formula1)
  e2 <- elementsInFormula(formula2)
  if (!identical(e1,e2)){
    e1 <- unique(append(e1,e2))
  }
  return(e1[order(e1, decreasing = decrease)])
}

#' @title Adding up two formulas, taking into account possible differing
#'  elements
#'
#' @param formula1 named numeric vector, example c(O = 2, C = 1)
#' @param formula2 named numeric vector, example c(H = 2, S = 1)
#'
#' @return a named numeric vector (formula)
#' @export
#'
#' @examples
#' addFormulas(waterFormula(), protonFormula())
#' addFormulas(waterFormula(), c(C=1, O=2))
addFormulas <- function(formula1, formula2){
  combinedNames <- elementsInFormulas(formula1, formula2)
  combined <- rep(0,length(combinedNames))
  names(combined) <- combinedNames
  unlist(lapply(names(combined),function(x){
    there <- which(names(formula1) == x)
    result <- 0
    if (!identical(there, integer(0))){
      result <- formula1[there]
    }
    there <- which(names(formula2) == x)
    if (!identical(there, integer(0))){
      result <- result + formula2[there]
    }
    return(result)
  }
  ))
}

#' @title custom operator for adding up formulas, to make calculating with
#'  formulas a little more clear
#' @param formula1 named numeric vector, example c(O = 2, C = 1)
#' @param formula2 named numeric vector, example c(H = 2, S = 1)
#' @examples
#' waterFormula() %f+% protonFormula()
#' waterFormula() %f+% c(C=1, O = 2)
#' c(H = 2, O = 1) %f+% c(S = 1, O = 2)
#'
#' @export
`%f+%` <- addFormulas


#' @title subtracting one formula from another, taking into account possible
#'  differing elements
#'
#' @param formula1 named numeric vector, example c(O = 2, C = 1); formula to be
#'  subtracted from
#' @param formula2 named numeric vector, example c(H = 2, S = 1); formula to
#'  subtract
#' @note There are no checks for negative values!
#' @return a named numeric vector (formula)
#' @export
#'
#' @examples
#' subtractFormulas(c(H = 2, O = 1), c(H = 1))
#' subtractFormulas(c(H = 2, O = 1), c(S = 1, O = 2))
subtractFormulas <- function(formula1, formula2){
  combinedNames <- elementsInFormulas(formula1, formula2)
  combined <- rep(0,length(combinedNames))
  names(combined) <- combinedNames
  unlist(lapply(names(combined),function(x){
    there <- which(names(formula1) == x)
    result <- 0
    if (!identical(there, integer(0))){
      result <- formula1[there]
    }
    there <- which(names(formula2) == x)
    if (!identical(there, integer(0))){
      result <- result - formula2[there]
    }
    return(result)
  }
  ))
}

#' @title custom operator for subtracting formulas from one another, to make
#'  calculating with formulas a little more clear
#' @param formula1 named numeric vector, example c(O = 2, C = 1); formula to be
#'  subtracted from
#' @param formula2 named numeric vector, example c(H = 2, S = 1); formula to
#'  subtract
#' @examples
#' c(H = 2, O = 1) %f-% c(H = 1)
#' c(H = 2, O = 1) %f-% c(S = 1, O = 2)
#'
#' @export
`%f-%` <- subtractFormulas

#' @title Add up a list of formulas
#'
#' @description Take a list of formulas and adds them all up
#'
#' @param formulas list of formulas
#'
#' @return a named numeric vector (formula)
#' @export
#'
#' @examples
#' addListFormulas(list(c(H = 2, O = 1),
#'                      c(H = 1),
#'                      c(H = 2, O = 1),
#'                      c(S = 1, O = 2)))
addListFormulas <- function(formulas){
  result <- formulas[[1]]
  for (i in 2:length(formulas)){
    result <- addFormulas(result, formulas[[i]])
  }
  return(result)
}

#' @title checks if formula is valid
#'
#' @note This is done via the check_chemform function from the package enviPat
#'
#' @param formula character vector or named numeric vector, representing the
#'  formula to be checked
#' @param string logical vector specifying if the formula is a character vector
#'  or not
#'
#' @return logical vector
#' @export
#'
#' @examples
#' glucose <- c(C=6, H=12, O=6)
#' validFormula(glucose)
#' formulaString(glucose)
#' validFormula(formulaString(glucose), string = TRUE)
validFormula <- function(formula, string = FALSE){
  if (!string){
    formula <- formulaString(formula = formula)
  }
  data(isotopes, envir = environment(), package = "enviPat")
  result <- enviPat::check_chemform(isotopes = isotopes, chemforms = formula)
  return(!result$warning)
}

#' @title calculates the neutral mono-isotopic mass of a formula
#'
#' @param formula named numeric vector, example c(O = 2, C = 1)
#' @param removeNA logical vector: what to do if any of the elements is NA. If
#'  TRUE, then remove before calculation, if FALSE, then do not remove
#' @param elementsInfo elements masses to be used, needs to be of class
#'  elements, default is elementsMonoisotopic(). The elementsAverage() function
#'  does not produce 100% correct answer for a lot of molecules due to the
#'  complex isotope patterns that emerge. In case average masses are needed it's
#'  better to use the enviPat option
#' @param enviPat logical argument that determines if the enviPat based
#'  calculations should be used. Default is FALSE. For monoisotopoc masses
#'  there is no difference, but for average masses of larger molecules (with
#'  complicated isotope patterns) it's highly recommended to use enviPat = TRUE
#'  with exact = FALSE
#' @param exact determines if the exact (TRUE, default) or the average (FALSE)
#'  mass is calculated (ignored if enviPat is FALSE)#'
#'
#' @return numeric vector
#' @export
#'
#' @examples
#' formulaToMass(c(H=2, O=1))
#' formulaToMass(c(H=2, O=1), elementsInfo = elementsAverage())
#' formulaToMass(c(C = 50, H=102))
#' formulaToMass(c(C = 50, H=102), elementsInfo = elementsAverage())
#' formulaToMass(c(C = 50, H=102), enviPat = TRUE)
#' formulaToMass(c(C = 50, H=102), enviPat = TRUE, exact = FALSE)
formulaToMass <- function(formula = NULL, removeNA = FALSE,
                          elementsInfo = elementsMonoisotopic(),
                          enviPat = FALSE, exact = TRUE){
  if (!is.null(formula)){
    if (sum(is.na(formula))>0){  # in case any of the elements is NA
      if (!removeNA){
        return(NA)
      } else {
        formula <- formula[!is.na(formula)]
      }
    }
    formula <- formula[formula > 0] # removeZeros
    if (length(formula) == 0){
      return(0)
    }
    if (!enviPat){
      # do not calculate via enviPat
      return(sum(unlist(lapply(names(formula),function(x){elementsInfo$getMass(x)})) * formula))
    } else {
      if (!validFormula(formula)){
        return(NA)
      } else {
        # calculate via enviPat
        data(isotopes, envir = environment(), package = "enviPat")
        if (!exact){
          suppressMessages(
            result <- enviPat::isopattern(isotopes = isotopes,
                                          chemforms = formula |> formulaString(),
                                          threshold = 0.01,
                                          plotit = FALSE,
                                          verbose = FALSE)[[1]]
          )
          return(stats::weighted.mean(result[,1], result[,2]))
        } else {
          result <- enviPat::check_chemform(isotopes = isotopes,
                                            chemforms = formula |> formulaString())$monoisotopic_mass
          if (result < 0){
            return(NA)
          } else {
            return(result)
          }
        }
      }
    }
  } else {
    return(NA)
  }
}

#' @title Calculates the m/z value of a charged/adducted ion
#'
#' @param mass numeric vector, (neutral) mass of the molecule
#' @param adducts numeric vector, number of adducts 'attached to' or 'removed
#'  from' the (originally neutral) molecule
#' @param adductFormula formula (named numeric vector) of the adduct
#' @param adductCharge numeric vector indicating the actual charge per adduct
#' @param elementsInfo elements masses to be used, needs to be of class
#'  elements, default is elementsMonoisotopic()
#'
#' @return numeric vector
#' @export
#'
#' @examples
#' # amino acid residue lysine + water
#' lysineMass <- formulaToMass(aminoAcidResidues()$getFormula("K") %f+% waterFormula())
#' lysineMass
#' # M+H+ : adducts =  1, adductFormula = protonFormula(), adductCharge = 1
#' # singly charged/protonated ion (ESI)
#' massToMz(lysineMass,
#'          adducts = 1,
#'          adductFormula = protonFormula(),
#'          adductCharge = 1)
#' # Doubly charged/protonated ion (ESI)
#' massToMz(lysineMass,
#'          adducts = 2,
#'          adductFormula = protonFormula(),
#'          adductCharge = 1)
#' # M-H- : single, negatively charged
#' massToMz(lysineMass,
#'          adducts = -1,
#'          adductFormula = protonFormula(),
#'          adductCharge = 1)
#' # M+ : singly positively charged (molecular) ion (EI)
#' massToMz(lysineMass,
#'          adducts = -1,
#'          adductFormula = electronFormula(),
#'          adductCharge = 1)
#' # M- : singly negatively charged (molecular) ion (EI)
#' massToMz(lysineMass,
#'          adducts = 1,
#'          adductFormula = electronFormula(),
#'          adductCharge = 1)
massToMz <- function(mass, adducts = 0, adductFormula = electronFormula(),
                     adductCharge = -1, elementsInfo = elementsMonoisotopic()){
  if (adducts == 0){
    return(mass)
  } else {
    massIon <- mass + (adducts * formulaToMass(adductFormula, elementsInfo = elementsInfo))
    # if charger is not an electron, then the (charge * mass(electron)) needs
    # to to be added/subtracted
    if (!identical(adductFormula,electronFormula())){
      massIon <- massIon - (abs(adductCharge) * adducts * formulaToMass(electronFormula(),
                                                                        elementsInfo = elementsInfo))
    }
    mz <- massIon / (abs(adducts * adductCharge))
    return(mz)
  }
}


#' @title Calculates the m/z value of a protonated ion (positive ESI)
#'
#' @description a wrapper around \code{\link{massToMz}} for positively charged,
#'  protonated ions in ESI
#'
#' @param mass numeric vector, (neutral) mass of the molecule
#' @param charge charge state
#' @param elementsInfo elements masses to be used, needs to be of class
#'  elements, default is elementsMonoisotopic()
#'
#' @return numeric vector
#' @export
#'
#' @examples
#' # amino acid residue lysine + water
#' lysineMass <- formulaToMass(aminoAcidResidues()$getFormula("K") %f+% waterFormula())
#' lysineMass
#' # M+H+ : adducts =  1, adductFormula = protonFormula(), adductCharge = 1
#' # singly charged/protonated ion (ESI)
#' massToMz(lysineMass,
#'          adducts = 1,
#'          adductFormula = protonFormula(),
#'          adductCharge = 1)
#' massToMzH(lysineMass)
massToMzH <- function(mass, charge = 1, elementsInfo = elementsMonoisotopic()){
  massToMz(mass, adducts = charge,
           adductFormula = protonFormula(), adductCharge = 1,
           elementsInfo = elementsInfo)
}


#' @title calculates the mass of the molecule in an ion
#'
#' @description essentially the reverse of the \code{\link{massToMz}}
#'
#' @param mz mass to charge ratio of the ion
#' @param adducts numeric vector, number of adducts 'attached to' or 'removed
#'  from' the (originally neutral) molecule
#' @param adductFormula formula (named numeric vector) of the adduct
#' @param adductCharge numeric vector indicating the actual charge per adduct
#' @param elementsInfo elements masses to be used, needs to be of class
#'  elements, default is elementsMonoisotopic()
#'
#' @return numeric vector
#' @export
#'
#' @examples
#' massToMz(mass = 174.1117, adductFormula = c(e=1), adducts = 2, adductCharge = -1) |>
#'  mzToMass(adductFormula = c(e=1), adducts = 2, adductCharge = -1)
#' massToMz(mass = 174.1117, adductFormula = c(H=1), adducts = 1, adductCharge = 1) |>
#'  mzToMass(adductFormula = c(H=1), adducts = 1, adductCharge = 1)
mzToMass <- function(mz, adducts = 0, adductFormula = electronFormula(), adductCharge = -1, elementsInfo = elementsMonoisotopic()){
  # deduct adducts
  result <- (mz * adducts * abs(adductCharge)) - (adducts * formulaToMass(adductFormula, elementsInfo = elementsInfo))
  # if charger is not electron
  if (!identical(adductFormula,electronFormula())){
    result <- result + (abs(adductCharge) * adducts * formulaToMass(electronFormula(),
                                                                    elementsInfo = elementsInfo))
  }
  return(result)
}

#' @title calculates the mass of the molecule in an ion (M+xH)x+
#'
#' @description a wrapper around \code{\link{mzToMass}} for positively charged,
#'  protonated ions in ESI
#'
#' @param mz numeric vector, mass to charge ratio of he ion
#' @param charge charge state
#' @param elementsInfo elements masses to be used, needs to be of class
#'  elements, default is elementsMonoisotopic()
#'
#' @return nuemric vector
#' @export
#'
#' @examples
#'  massToMzH(mass = 174.1117, charge = 2) |> mzHToMass(charge = 2)
#'  massToMzH(mass = 174.1117, charge = 1) |> mzHToMass(charge = 1)
mzHToMass <- function(mz, charge = 1, elementsInfo = elementsMonoisotopic()){
  mzToMass(mz = mz, adducts = charge,
           adductFormula = protonFormula(), adductCharge = 1,
           elementsInfo = elementsInfo)
}

# ---- translations from and to different formats ----

#' @title  Translates regular formula format into a character vector, eg
#'  C6H12O6
#'
#' @param formula named numeric vector, example c(O = 2, C = 1)
#' @param removeSingle if TRUE then for elements that are present in the
#'  formula only a single time, the number (1) will not be included. Default is
#'  FALSE. See also examples
#' @param useMarkdown default = FALSE. If TRUE, the it will use HTML/Markdown
#'  codes <sub> in the formulas, which can be used with the library 'gt' to
#'  generate 'proper' notation for chemical formulas (numbers in subscript)
#'
#' @return character vector
#' @export
#'
#' @examples
#' formulaString(c(C=6, H=12, O=6))
#' formulaString(c(H=3,O=4,P=1))
#' formulaString(c(H=3,O=4,P=1), removeSingle = TRUE)
#' formulaString(c(H=2, O=1))
#' formulaString(c(H=2, O=1), removeSingle = TRUE)
formulaString <- function(formula, removeSingle = FALSE, useMarkdown = FALSE){
  for (counter in 1:length(formula)){
    if (stringr::str_count(names(formula)[counter],
                           pattern = "\\d+[:alpha:]+")>0){
      # eg 13C
      locations <- stringr::str_locate(names(formula)[counter],
                                       pattern = "\\d+[:alpha:]+")
      names(formula)[counter] <- paste(
        c("[",
          substr(names(formula)[counter],
                 locations[1],
                 locations[2]-1),
          "]",
          substr(names(formula)[counter],
                 locations[2], nchar(names(formula)[counter]))
        ), collapse = "")
    }
  }
  result <- paste(
    unlist(
      lapply(1:length(formula),
             function(x){paste(names(formula)[x],
                               ifelse(removeSingle & (formula[x] == 1),
                                      "",
                                      paste(c(ifelse(useMarkdown,"<sub>",""),
                                              toString(formula[x]),
                                              ifelse(useMarkdown,"</sub>","")),
                                            collapse = "")),
                               sep = "")})),
    collapse ="")

  if (useMarkdown){
    result <- gsub(result, pattern = "\\[", replacement = "<sup>")
    result <- gsub(result, pattern = "\\]", replacement = "</sup>")
  }
  return(result)
}

#' @title  Translates a character vector formula, eg 'C6H12O6' to a regular
#'  formula c(C=6, H=12, O=6)
#'
#' @note it's imperative that every element has a number (count), otherwise this
#'  function is highly likely to malfunction and return NA
#'
#' @param string character vector, format eg: 'C6H12O6'
#'
#' @return formula of format c(H=2, O=1)
#' @export
#'
#' @examples
#' stringFormula("H3O4P1")
#' stringFormula("C6H12O6")
stringFormula <- function(string){
  parts <- stringr::str_extract_all(string,
                                    pattern = "((\\[\\d+\\]){0,1}[:alpha:]+\\-{0,1}\\d*)")[[1]]
  elementsList <- as.character()
  elementsCount <- as.numeric()
  for (counter in 1:length(parts)){
    theElement <- stringr::str_extract_all(parts[counter],
                                           pattern= "(\\[\\d+\\]){0,1}[:alpha:]+")[[1]]
    elementCounted <- stringr::str_extract_all(parts[counter],
                                               pattern= "\\-{0,1}\\d+")[[1]]
    if (length(elementCounted) > 0){
      # in case of isotope eg 13C5
      elementCounted <- as.numeric(elementCounted[length(elementCounted)])
    } else {
      # note: this is highly likely to malfunction, ONLY works when last element
      # in the list does not have a number
      elementCounted <- 1
    }
    theElement <- stringr::str_replace_all(theElement, pattern = "\\[", replacement = "")
    theElement <- stringr::str_replace_all(theElement, pattern = "\\]", replacement = "")
    elementsList <- append(elementsList, theElement)
    elementsCount <- append(elementsCount, elementCounted)
  }
  names(elementsCount) <- elementsList
  return(elementsCount)
}

#' @title Translates a character vector formula, eg 'C6H12O6' to a regular
#'  formula c(C=6, H=12, O=6)
#'
#' @note this function is an improved version of stringFormula(). Now every elements
#'  with count 1 can have the number omitted. However, the function depends on 'correct'
#'  elements (first letter is uppercase, second letter is lowercase). This function also
#'  allows for the presence of isotopes, eg '[13]C' or '[2]H2O'
#'
#' @param string character vector, format eg: 'C6H12O6'
#'
#' @return formula of format c(H=2, O=1)
#' @export
#'
#' @examples
#' stringToFormula("H3O4P1")
#' stringToFormula("C6H12O6")
#' stringToFormula("C6H5Br")
#' stringToFormula("[13]C6H12O5[18]O")
stringToFormula <- function(string){
  if (nchar(string)==0){
    return(emptyFormula())
  }
  if (length(string)==1){
    findElements <- stringr::str_locate_all(string,
                                            pattern = "(\\[{0,1}[:digit:]+\\]{0,1}){0,1}[:upper:]{1}[:lower:]{0,1}\\d*")[[1]]
    result <- as.numeric()
    for (counter in 1:nrow(findElements)){
      tempStr <- stringr::str_sub(string,
                                  start = findElements[counter, "start"],
                                  end = findElements[counter, "end"])
      result[counter] <- as.numeric(stringr::str_extract(tempStr, pattern = "\\d+$"))
      if (is.na(result[counter])){
        result[counter] <- 1
      }
      names(result)[counter] <- stringr::str_extract(tempStr, pattern = ".*[:alpha:]+")
    }
    return(result)
  }
  return(lapply(string, stringToFormula))
}

#' @title translates an cdkFormula object to a 'regular' formula format
#'
#' @param cdkformula an object of type rcdkFormula
#'
#' @note This function does not deal with the charge state which is possibly
#'  defined in the rcdkFormula object
#'
#' @return formula of format c(H=2, O=1)
#' @export
#'
#' @examples
#' glucose <- rcdk::get.formula("C6H12O6")
#' rcdkFormula(glucose)
#' glucoseAdductIon <- rcdk::get.formula("C6H12O6Na1", charge = 1)
#' glucoseAdductIon
#' # to get to the same m/z value
#' glucoseAdductIon |> massSpectrometryR::rcdkFormula() |> formulaToMass() |> massToMz(adducts = -1)
rcdkFormula <- function(cdkformula){
  tempF <- as.data.frame(cdkformula@isotopes)[,c("isoto", "number")]
  result <- as.numeric(tempF$number)
  names(result) <- tempF$isoto
  return(result)
}

#' @title translates a proteome Discoverer (Thermo Scientific) elements formula string
#'  to a formula as used by this package
#'
#' @param pdFormula a character vector. Formula in a format as used by proteome
#'  discoverer software
#'
#' @return formula of format c(H=2, O=1)
#' @export
#'
#' @examples
#' glucose <- pdToFormula("C(6) H(12) O(6)")
#' glucose
#' water <- pdToFormula("H(2) O")
#' water
pdToFormula <- function(pdFormula){
  if (length(pdFormula) == 1){
    pdFormula <- stringr::str_split(pdFormula, pattern = " ")[[1]]
    findWOnumber <- !grepl(pdFormula, pattern = "\\(\\d+\\)")
    pdFormula[findWOnumber] <- paste(pdFormula[findWOnumber],"(1)", sep = "")
    pdFormula <- stringr::str_split(
      stringr::str_replace(pdFormula,
                           pattern = "\\)",
                           replacement = ""),
      pattern = "\\(")
    pdFormula <- unlist(pdFormula)
    result <- as.integer(pdFormula[seq(2, length(pdFormula),2)])
    names(result) <- pdFormula[seq(1, length(pdFormula),2)]
  } else {
    result <- purrr::map(pdFormula, ~pdTomsFormula(.x))
    names(result) <- names(pdFormula)
  }
  return(result)
}

# ---- Standard Formulas ----

#' @title generates a pre-defined formula for proton
#' @return a named numeric vector (formula)
#' @examples
#' print(protonFormula())
#' @export
protonFormula <- function(){
  return(c(H = 1))
}

#' @title generates a pre-defined formula for water
#' @return a named numeric vector (formula)
#' @examples
#' print(waterFormula())
#' @export
waterFormula <- function(){
  return(c(H = 2, O = 1))
}

#' @title generates a pre-defined formula for electron
#' @note this is used for calculations
#' @return a named numeric vector (formula)
#' @examples
#' print(electronFormula())
#' @export
electronFormula <- function(){
  return(c(e = 1))
}

#' @title generates an empty pre-defined formula
#' @note this can be used for calculations and setup of unknowns
#' @return a named numeric vector (formula)
#' @examples
#' print(emptyFormula())
#' @export
emptyFormula <- function(){
  return(c(C = 0, H =  0, N = 0, O = 0, S = 0))
}

# ---- chemicals ----

#' R6 Class representing a set of chemicals
#'
#' Note: this class is meant to be used for classes of compounds, eg amino acids
#'
#' Also: his class is  meant as a base class to be expanded via inheritance
#'
#' Warning: all chemicals inside this object should be unique (names, letters &
#'  shorts)
#'
#' @description
#' Every chemical inside the object has a name, letter, short and a formula.
#'  The first 3 can be any length of string (though the letter and short field
#'   should be maximum length (nchar) 1 and 2-4 respectively). Formula should be
#'   in the form of a named numeric with the names representing elements and the
#'   values themselves being the number of atoms of that element,
#'   eg c(C = 3, H =  5, N = 1, O = 1, S = 0)
#'
#' @examples
#' estrogens <- chemicals$new(letters = c("1","2","3","4"),
#'                            shorts = c("E1","E2","E3","E4"),
#'                            names = c("Estrone","Estradiol",
#'                                      "Estriol","Estetrol"),
#'                            formulas = list(c(C=18, H=22, O=2),
#'                                            c(C=18, H=24, O=2),
#'                                            c(C=18, H=24, O=3),
#'                                            c(C=18, H=24, O=4)))
#'
#' @export
chemicals <- R6::R6Class("chemicals",
                      private = list(
                        # letters can be simply a number of some sorts, but
                        # actaully meant to store eg aminoa acid letters,
                        # eg G for Glycine
                        letters_ = NA,
                        # shorts the short names for the chemicals,
                        # eg Ala for Alanine
                        shorts_ = NA,
                        # names the actual names of the chemicals
                        names_ = NA,
                        # formulas for the chemicals,
                        # eg c(C = 6, H = 12, N = 4, O = 1, S = 0) for Arginine
                        formulas_ = NA
                      ),
                      public = list(
                        #' @description
                        #' Create a new chemicals object
                        #' @param letters character vector specifying the letters (or numbers or whatever) for the chemicals.
                        #'  In case of amino acids it should be eg "A" for Alanine, "G" for Glycine, etc etc
                        #' @param shorts  character vector specifying the short names for the chemicals, eg Ala for Alanine
                        #' @param names character vector specifying the names of the chemicals
                        #' @param formulas list of named numeric vectors specifying the formulas of the chemicals,
                        #'  eg c(C = 6, H = 12, N = 4, O = 1, S = 0) for Arginine
                        #' @return a new 'chemical' object
                        initialize = function(letters, shorts,
                                              names, formulas){
                        stopifnot(length(letters) == length(shorts) &
                                    length(letters) == length(names) &
                                    length(letters) == length(formulas))
                          private$letters_  <- letters
                          private$names_    <- names
                          private$shorts_   <- shorts
                          private$formulas_ <- formulas
                          invisible(self)
                        },
                        #' @description
                        #' For printing purposes: prints a table of the chemicals with columns letter, name & short
                        #' @param ... no arguments, the function takes care of printing
                        print = function(...){
                          print(dplyr::tibble(self$table))
                        },
                        #' @description Retrieves the formula of one of the compounds in the object
                        #' @param which1 specifies which chemical should be retrieved. Which the number (row number in the chemicals
                        #'  table), or the name, letter or short as a character vector. The way this is set up, it doesn't matter
                        #'  whether capital or non-capital letters are used, since all is converted to upper case before comparing
                        #'  with what's in the chemical table
                        #' @return a formula in the shape of a named numeric vector, eg c(C = 6, H = 12, N = 4, O = 1, S = 0)
                        getFormula = function(which1){
                          stopifnot(!is.null(which1) & !is.na(which1))
                          # number
                          if (which1 %in% 1:self$number){
                            return(unlist(private$formulas_[which1]))
                          } else {
                            which1 <- toupper(which1)
                            # letter
                            temporary <- which(toupper(private$letters_) == which1)
                            if (!(identical(temporary, integer(0)))){
                              return(unlist(private$formulas_[temporary]))
                            } else {
                              # 3-letter short
                              temporary <- which(toupper(private$shorts_) == which1)
                              if (!(identical(temporary, integer(0)))){
                                return(unlist(private$formulas_[temporary]))
                              } else { # full name
                                temporary <- which(toupper(private$names_) == which1)
                                if (!(identical(temporary, integer(0)))){
                                  return(unlist(private$formulas_[temporary]))
                                } else {
                                  return(NA)
                                }
                              }
                            }
                          }
                        }
                      ),
                      active = list(
                        #' @field number retrieve the number of compounds present in the object, read only
                        number = function(value){
                          if (missing(value)){
                            return(length(private$formulas_)) # assumption is that all 4 fields are same length
                          } else {
                            # nothing, read only!
                          }
                        },
                        #' @field letters to access the letters of the compounds in the object
                        letters = function(value){
                          if (missing(value)){
                            return(private$letters_)
                          } else {
                            private$letters_ <- value
                          }
                        },
                        #' @field names to access the names of the compounds in the object
                        names = function(value){
                          if (missing(value)){
                            return(private$names_)
                          } else {
                            private$names_ <- value
                          }
                        },
                        #' @field shorts to access the shorts of the compounds in the object
                        shorts = function(value){
                          if (missing(value)){
                            return(private$shorts_)
                          } else {
                            private$shorts_ <- value
                          }
                        },
                        #' @field formulas to access the formulas of the compounds in the object
                        formulas = function(value){
                          if (missing(value)){
                            return(private$formulas_)
                          } else {
                            private$formulas_ <- value
                          }
                        },
                        #' @field table retrieves all info on the compounds in data.frame format, read only
                        table = function(value){
                          if (missing(value)){
                            df_ <- data.frame(
                              letter = private$letters_,
                              name = private$names_,
                              short = private$shorts_)
                            df_$formula <- private$formulas_
                            return(df_)
                          } else {
                            # nothing, read only
                          }
                        }
                      )
)


#' R6 Class representing a set of amino acids
#'
#' Note: this class is meant to be used only for amino acids and such
#'
#' @description
#' R6 Class representing a set of amino acids. It adds three functions to quickly
#'  switch between different writing 'styles' of peptides
#'
#' @examples
#' aminoAcidResidues()$getShort("L")
#' aminoAcidResidues()$getShort("Leu")
#' aminoAcidResidues()$getName("L")
#' aminoAcidResidues()$getName("Leu")
#' aminoAcidResidues()$translatePeptide("Asp-Arg-Val-Tyr-Ile-His-Pro-Phe-His-Leu",
#'  from1to3 = TRUE, splitCharacter ="-")
#' aminoAcidResidues()$translatePeptide("DRVYIHPFHL", joinCharacter = "-")
#'
#' @export
aminoAcidClass <- R6::R6Class("aminoacids",
                              inherit = chemicals,
                              public = (
                                list(
                                #' @description
                                #' Function to retrieve the full name of an amino acid via the letter or shorts
                                #' @param searchString either a 1- or 3- letter character vector
                                #' @param checkCase default = TRUE. If false, the function will ignore the case
                                #'  the searchString argument
                                #'
                                #' @return character vector, name of the aminoacid
                                getName = function(searchString, checkCase = TRUE){
                                  if (checkCase){
                                    if (nchar(searchString) == 1){
                                      if (searchString %in% private$letters_){
                                        return(private$names_[which(private$letters_ == searchString)])
                                      }
                                    } else {
                                      if (nchar(searchString) == 3){
                                        if (searchString %in% private$shorts_){
                                          return(private$names_[which(private$shorts_ == searchString)])
                                        }
                                      }
                                    }
                                  } else {
                                    if (nchar(searchString) == 1){
                                      if (toupper(searchString) %in% toupper(private$letters_)){
                                        return(private$names_[which(toupper(private$letters_) == toupper(searchString))])
                                      }
                                    } else {
                                      if (nchar(searchString) == 3){
                                        if (toupper(searchString) %in% toupper(private$shorts_)){
                                          return(private$names_[which(toupper(private$shorts_) == toupper(searchString))])
                                        }
                                      }
                                    }
                                  }
                                  return(NA)
                                },
                                #' @description
                                #' Function to retrieve either the 1- or 3- letter code of an amino acid
                                #' @param searchString either a 1- or 3- letter character vector: if 1-letter
                                #'  than the corresponding 3-letter character vector will be returned and vice
                                #'  versa
                                #' @param checkCase default = TRUE. If false, the function will ignore the case
                                #'  the searchString argument
                                #'
                                #' @return character vector, 1- or 3- letter code of the aminoacid
                                getShort = function(searchString, checkCase = TRUE){
                                  if (checkCase){
                                    if (nchar(searchString) == 1){
                                      if (searchString %in% private$letters_){
                                        return(private$shorts_[which(private$letters_ == searchString)])
                                      }
                                    } else {
                                      if (nchar(searchString) == 3){
                                        if (searchString %in% private$shorts_){
                                          return(private$letters_[which(private$shorts_ == searchString)])
                                        }
                                      }
                                    }
                                  } else {
                                    if (nchar(searchString) == 1){
                                      if (toupper(searchString) %in% toupper(private$letters_)){
                                        return(private$shorts_[which(toupper(private$letters_) == toupper(searchString))])
                                      }
                                    } else {
                                      if (nchar(searchString) == 3){
                                        if (toupper(searchString) %in% toupper(private$shorts_)){
                                          return(private$letters_[which(toupper(private$shorts_) == toupper(searchString))])
                                        }
                                      }
                                    }
                                  }
                                  return(NA)
                                },
                                #' @description
                                #' Translates a amino acid sequence from 1-letter codes to 3-letter codes and
                                #'  vice versa
                                #'
                                #' @param sequence character vector: amino acid sequence in 1-letter or
                                #'  3-letter codes
                                #' @param from1to3 logical vector: if TRUE, then translation will be from
                                #'  1-letter code to 3=letter code. If FALSE, then vice versa. Default = FALSE
                                #' @param splitCharacter character vector specifying the character(s) between
                                #'  the 1- or 3-letter codes in the sequence. Default NA (same as "")
                                #' @param joinCharacter character vector specifying the character(s) between
                                #'  the translated codes. Default NA (same as "")
                                #' @param checkCase default = TRUE. If false, the function will ignore the case
                                #'  the sequence
                                #'
                                #' @return character vector, sequence in either 1- or 3-letter codes
                                translatePeptide = function(sequence,
                                                            from1to3 = FALSE, splitCharacter = NA, joinCharacter = NA,
                                                            checkCase = TRUE){
                                  if (identical(splitCharacter, NA)){
                                    splitCharacter =""
                                  }
                                  if (from1to3){
                                    sequence <- stringr::str_split(sequence, pattern = splitCharacter)[[1]]
                                  } else {
                                    if (splitCharacter == ""){
                                      sequence <- stringr::str_extract_all(sequence, pattern = "[:alpha:]{3}")[[1]]
                                    } else {
                                      sequence <- stringr::str_split(sequence, pattern = splitCharacter)[[1]]
                                    }
                                  }
                                  return(paste(purrr::map_chr(sequence,
                                                              ~self$getShort(.x,
                                                                             checkCase = checkCase)),
                                               collapse = ifelse(identical(joinCharacter, NA),
                                                                 "",
                                                                 joinCharacter)))
                                }
                               )
                              )
)

#' @title Generates a pre-defined object which contains info on 'normal' amino
#'  acid residues
#' @note The formulas in the object are amino acid residues as they are present in proteins.
#'  To get the actual formula of the amino acid in its 'free' form, add c(H=2, O=1) (water)
#' @note this object is used in all protein calculations in this package
#'
#' @return a R6 object of class 'chemicals'
#' @examples
#' print(aminoAcidResidues())
#'
#' @export
aminoAcidResidues <- function(){
  aminoAcidClass$new(letters =  c("A","R","N","D",
                                  "E","Q","G","H",
                                  "I","L","K","M",
                                  "F","P","S","T",
                                  "W","Y","V","C"),
                     shorts   = c("Ala","Arg","Asn","Asp",
                                  "Glu","Gln","Gly","His",
                                  "Ile","Leu","Lys","Met",
                                  "Phe","Pro","Ser","Thr",
                                  "Trp","Tyr","Val","Cys"),
                   names    = c("Alanine","Arginine","Asparagine", "Aspartate",
                                "Glutamate","Glutamine","Glycine","Histidine",
                                "Isoleucine","Leucine","Lysine","Methionine",
                                "Phenylalanine", "Proline","Serine","Threonine",
                                "Tryptophane", "Tyrosine","Valine", "Cysteine"),
                     formulas = list(c(C = 3, H =  5, N = 1, O = 1, S = 0),
                                     c(C = 6, H = 12, N = 4, O = 1, S = 0),
                                     c(C = 4, H =  6, N = 2, O = 2, S = 0),
                                     c(C = 4, H =  5, N = 1, O = 3, S = 0),
                                     c(C = 5, H =  7, N = 1, O = 3, S = 0),
                                     c(C = 5, H =  8, N = 2, O = 2, S = 0),
                                     c(C = 2, H =  3, N = 1, O = 1, S = 0),
                                     c(C = 6, H =  7, N = 3, O = 1, S = 0),
                                     c(C = 6, H = 11, N = 1, O = 1, S = 0),
                                     c(C = 6, H = 11, N = 1, O = 1, S = 0),
                                     c(C = 6, H = 12, N = 2, O = 1, S = 0),
                                     c(C = 5, H =  9, N = 1, O = 1, S = 1),
                                     c(C = 9, H =  9, N = 1, O = 1, S = 0),
                                     c(C = 5, H =  7, N = 1, O = 1, S = 0),
                                     c(C = 3, H =  5, N = 1, O = 2, S = 0),
                                     c(C = 4, H =  7, N = 1, O = 2, S = 0),
                                     c(C = 11,H = 10, N = 2, O = 1, S = 0),
                                     c(C = 9, H =  9, N = 1, O = 2, S = 0),
                                     c(C = 5, H =  9, N = 1, O = 1, S = 0),
                                     c(C = 3, H =  5, N = 1, O = 1, S = 1)))
}
