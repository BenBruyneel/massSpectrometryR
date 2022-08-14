# ---- General ----

#' helper function to determine if object == whichClass or a descendant
#'  of whichClass.
#'
#' @param object a data object of some class
#' @param whichClass character string: class name to be tested
#' @return TRUE or FALSE
#' @noRd
is.Class <- function(object, whichClass){
  return(whichClass %in% class(object))
}

#' @title Generates a pre-defined (incomplete) table of names of fragments for
#'  the ions resulting when fragmenting a peptide in MS
#' @note will possibly be removed in the future
#' 
#' @return a data.frame of two columns: 'name' and 'series name' of fragments
#' @export
peptideFragments <- function(){
  return(
    data.frame(name = c('aions', 'aionsH2O', 'aionsNH3',
                        'bions', 'bionsH2O', 'bionsNH3',
                        'bionsH2Oadd', 'cions', 'xions',
                        'yions', 'yionsH2O', 'yionsNH3',
                        'zions'),
               series = c('a', 'a-H2O', 'a-NH3',
                          'b', 'b-H2O', 'b-NH3',
                          'b+H2O', 'c', 'x',
                          'y', 'y-H2O', 'y-NH3',
                          'z'))
  )
}


# ---- general - calculations ----

#' @title peptideToFormula
#' @description gives formula of a peptide string c(C=x,H=y,etc etc)
#' 
#' @note does not check for non-amino acid letters, modifications cannot be
#'  specified
#'
#' @param peptide character vector specifying the sequence of amino acids in a
#'  peptide
#' @param aminoAcids R6 object of type 'chemicals' with the amino acid info,
#'  default = aminoAcidResidues()
#'
#' @return a numeric vector
#' @export
#'
#' @examples
#' peptideToFormula("SAMPLER")
peptideToFormula <- function(peptide, aminoAcids = aminoAcidResidues()){
  peptide <- unlist(stringr::str_split(toupper(peptide),pattern = ""))
  temporary <- lapply(peptide, aminoAcids$getFormula)
  # need to add 1 x water
  summ <- waterFormula()
  for (cntr in 1:length(temporary)){
    summ <- addFormulas(summ,temporary[[cntr]])
  }
  return(summ)
}


#' @title peptideMassToMz
#' @description gives ion m/z of a peptide string c(C=x,H=y,etc etc)
#'
#' @note does not check for non-amino acid letters, modifications cannot be
#'  specified
#'  
#' @param peptide character vector specifying the sequence of amino acids in a
#'  peptide
#' @param charge numeric vector specifying the charge of the peptide ion
#' @param aminoAcids R6 object of type 'chemicals' with the amino acid info,
#'  default = aminoAcidResidues()
#' @param elementsInfo R6 object of type 'elements' with the elements masses
#'  info, default = elementsMonoisotopic()
#'
#' @return a numeric vector
#' @export
#'
#' @examples
#' peptideMassToMz("SAMPLER")
#' peptideMassToMz("SAMPLER", charge = 2)
#' peptideMassToMz("SAMPLER", elementsInfo = elementsAverage())
peptideMassToMz <- function(peptide, charge = 1,
                            aminoAcids = aminoAcidResidues(),
                            elementsInfo = elementsMonoisotopic()){
  peptideToFormula(peptide, aminoAcids = aminoAcids) |>
    formulaToMass(elementsInfo = elementsInfo) |>
    massToMzH(charge = charge, elementsInfo = elementsInfo)
}


#' @title counts the occurence of a amino acid (sequence) in another amino acid
#'  sequence
#'
#' @param thePeptide character vector, the peptide to be searched 
#' @param searchPeptide character vector, the amino acid sequence to search for
#' @param doNotSplice if FALSE the all characters in the searchPeptide are
#'  searched individually.If TRUE then the searchPeptide is searched as a whole.
#'  Default = TRUE
#' @param upper convert both thePeptide & searchPeptides to uppercase before
#'  searching
#'
#' @return numeric vector
#' @export
#'
#' @examples
#' peptideCount("SAMPLER", "P")
#' peptideCount("SAMPLER", "PLER", doNotSplice = TRUE)
#' peptideCount("SAMPLER", "PLER", doNotSplice = FALSE)
peptideCount <- function(thePeptide = NA, searchPeptide = NA, doNotSplice = TRUE, upper = TRUE){
  if (!(is.na(thePeptide) | (identical(searchPeptide,NA)))){
    if (upper){
      thePeptide <- toupper(thePeptide)
      searchPeptide <- toupper(searchPeptide)        
    }
    if (length(searchPeptide) == 1){
      if (doNotSplice){
        return(stringr::str_count(thePeptide, searchPeptide))
      } else {
        return(unlist(lapply(unlist(stringr::str_split(searchPeptide,"")),function(x){peptideCount(thePeptide,x)})))
      }
    } else {
      return(unlist(lapply(searchPeptide, function(x){peptideCount(thePeptide,x)})))
    }
  } else {
    return(NA)
  }
}

# ---- modifications ----

#' R6 Class representing a set of modifications for the aminoacids in peptides
#' 
#' @description 
#' Every modification inside the object has a name, position, fixed (flag),
#'  gain (formula), loss (formula) and category.
#'  
#' 'name' is a character vector.
#' 
#'  Position is a character vector specifying the amino acids which always have
#'   the modification (fixed = TRUE) or can have the modification (fixed =
#'   FALSE). More than one amino acid can be specified, eg NQ (for  Asparagine &
#'   glutamine). Please note that currently the package does NOT support
#'   'exotic' amino acids (eg Selenocysteine) or 'combination' letters, such as
#'   'J' (Leucine or Isoleucine). For the C- and N-terminus, use 'C_Term' or
#'   'N_Term' for position.
#'   
#'  Gain and loss specify what is lost and/or gained when a amino acid is
#'   modified. For example Carbamidomethylation of Cysteine has both a loss
#'   formula c(H=1) and a gain formula c(C=2, H=4, N=1, O=1); obviously this
#'   could also be defined as: loss formula = emptyFormula(), gain formula =
#'   c(C=2, H=3, N=1, O=1).
#'   
#'   For the category field (character vector) there is no real 'rule' on how to
#'   classify modifications. I usually stick to the categorisation of Mascot or
#'   Sequest.
#'   
#' @note the logical vector fixed is very important. If TRUE, then a
#'  modification is considered to be always present, if FALSE then its presence
#'  is optional. 
#'   
#' @examples 
#' aaModifications <- modifications$new()
#' aaModifications$addTable(
#'   tibble::tibble(
#'     name = c("Carbamidomethyl (C)",
#'              "Carboxymethyl (C)",
#'              "Oxidation (M)"),
#'     position = c("C","C","M"),
#'     fixed = c(TRUE,TRUE,FALSE),
#'     gain = list(c(C = 2, H = 4, N = 1, O = 1),
#'                 c(C = 2, H = 3, N = 0, O = 2),
#'                 c(C = 0, H = 0, N = 0, O = 1)),
#'     loss = list(c(protonFormula()),
#'                 c(protonFormula()),
#'                 c(emptyFormula())),
#'     category = c("Cys-state","Cys-state","Preparation Artefact")
#'   )
#' )
#' aaModifications
#' aaModifications$add(name = "Deamidation",
#'                     position = "NQ",
#'                     fixed = FALSE,
#'                     gain = c(O = 1, H = 1),
#'                     loss = c(N = 1, H =2),
#'                     category = "Preparation Artefact")
#' aaModifications
#' @export
modifications <- R6Class("modifications",
                         private = list(
                           # table has 6 columns or is completely empty, see
                           # description & examples
                           modificationTable = tibble::tibble()
                         ),
                         public = list(
                           #' @description 
                           #' Create a new modifications object
                           #' @param data default = NA. If not NA, then should be a tibble with 6 columns: name, position,
                           #'  fixed, gain, loss and category. This is checked, but the contents of each column are
                           #'  not checked.
                           initialize = function(data = NA){
                             self$addTable(data = data)
                             invisible(self)
                           },
                           #' @description
                           #' For printing purposes: prints a table of the modifications
                           #' @param ... no arguments, the function takes care of printing
                           print = function(...){
                             if (identical(private$modificationTable,
                                           tibble::tibble())){
                               print(NA)
                             } else {
                               if (sum(private$modificationTable$fixed)>0){
                                 cat("\n")
                                 cat("Fixed Modifications\n")
                                 print(private$modificationTable[private$modificationTable$fixed,])
                               }
                               if (sum(!private$modificationTable$fixed)>0){
                                 cat("\n")
                                 cat("Variable Modifications\n")
                                 print(private$modificationTable[!private$modificationTable$fixed,])
                               }
                             }
                           },
                           #' @description 
                           #'  Adds a single modification
                           #'
                           #' @param name character vector
                           #' @param position character vector, should be a valid amino acid residue
                           #' @param fixed logical vector, specifies whether the modification is fixed (TRUE) or
                           #'  dynamic (FALSE)
                           #' @param gain named numeric vector (formula) that specifies what (atoms) are gained
                           #'  when a modification is applied to an amino acid
                           #' @param loss named numeric vector (formula) that specifies what (atoms) are lost
                           #'  when a modification is applied to an amino acid 
                           #' @param category character vector. Not rigidly defined: for user to be able to
                           #'  select/filter etc which type of modifications to use
                           #'
                           #' @return nothing
                           #' @export
                           add = function(name, position, fixed,
                                          gain, loss,
                                          category){
                             private$modificationTable <- dplyr::bind_rows(
                               private$modificationTable,
                               tibble::tibble(name = name,
                                              position = position,
                                              fixed = fixed,
                                              gain = list(gain),
                                              loss = list(loss),
                                              category = category))
                             invisible(self)
                           },
                           #' @description 
                           #'  Add (a set of) modifications via a tibble
                           #'
                           #' @param data default = NA. If not NA, then should be a tibble with 6 columns: name, position,
                           #'  fixed, gain, loss and category. This is checked, but the contents of each column are
                           #'  not checked.
                           #'
                           #' @return nothing
                           #' @export
                           addTable = function(data = NA){
                             if (!identical(data, NA)){
                               if (is.Class(data, "tbl_df")){  # since tibble itself is not a class name apparently
                                 if (sum(unlist(lapply(colnames(data),
                                                       function(x) {
                                                         !( x %in% c("name",
                                                                     "position",
                                                                     "fixed",
                                                                     "gain",
                                                                     "loss",
                                                                     "category"
                                                         ))}))) == 0){
                                  private$modificationTable <- dplyr::bind_rows(
                                    private$modificationTable,
                                    data)
                                 }
                               }
                             }
                             invisible(self)
                           }
                         ),
                         active = list(
                           #' @field number retrieve the number of modifications present in the object, read only
                           number = function(value){
                             if (missing(value)){
                               return(nrow(private$modificationTable))
                             } else {
                               # do nothing, read only
                             }
                           },
                           #' @field fixed retrieve a table of the fixed modifications in the modification table, read only
                           fixed = function(value){
                             if (missing(value)){
                               return(private$modificationTable[private$modificationTable$fixed,])
                             } else {
                               # do nothing, read only
                             }
                           },
                           #' @field variable retrieve a table of the variable modifications in the modification table, read only
                           variable = function(value){
                             if (missing(value)){
                               return(private$modificationTable[!private$modificationTable$fixed,])
                             } else {
                               # do nothing, read only
                             }
                           },
                           #' @field table to access the table of modifications directly
                           table = function(value){
                             if (missing(value)){
                               return(private$modificationTable)
                             } else {
                               private$modificationTable <- value
                             }
                           }
                         )
)                        

#' @title Returns a pre-defined object which contains info on some common
#'  amino acid modifications
#'  
#' @note the resulting modification table cannot be used immediately: there is
#'  two times a fixed modification for Cysteine amino acids. Remove one of them
#'  to prevent errors when using peptide calculations
#'  
#' @return An object of class modifications containing info on amino acid
#'  modifications
#' @examples
#' print(aminoAcidModifications)
#' 
#' @export
aminoAcidModifications <- function(){
  return(
    modifications$new(
      data = tibble::tibble(
        name = c("Carbamidomethyl (C)",
                 "Carboxymethyl (C)",
                 "Oxidation (M)",
                 "Deamidated (NQ)",
                 "Oxidation (W)",
                 "Oxidation Double (W)"),
        position = c("C","C","M","NQ","W","W"),
        fixed = c(TRUE,TRUE,FALSE,FALSE,FALSE,FALSE),
        gain = list(c(C = 2, H = 4, N = 1, O = 1),
                    c(C = 2, H = 3, N = 0, O = 2),
                    c(C = 0, H = 0, N = 0, O = 1),
                    c(O = 1, H = 1),
                    c(O = 1),
                    c(O = 2)),
        loss = list(protonFormula(),
                    protonFormula(),
                    emptyFormula(),
                    c(N = 1, H =2),
                    emptyFormula(),
                    emptyFormula()),
        category = c("Cys-state",
                     "Cys-state",
                     "Preparation Artefact",
                     "Preparation Artefact",
                     "Preparation Artefact",
                     "Preparation Artefact")
      )
    )
  )
}

# ---- peptide ----

#' R6 Class representing a (single) peptide
#' 
#' @description
#'  Contains two character vectors: one representing the amino acid sequence,
#'  and a second conatining info on the positions of 'variable' modifications.
#'  The object also contains a modification table specifying the 'fixed' amd
#'  'variable' modifications.
#' @examples
#' testPeptide <- peptide$new(sequence = "SAMPLER",
#'                            modificationTable = aminoAcidModifications()$table,
#'                            variableModifications = "0010000")
#' testPeptide
#' testPeptide$formula()                          
#' @export  
peptide <- R6Class("peptide",
                   private = list(
                     peptideString = NA,
                     # variable modifications
                     # format = string of same length as peptideString, similair to how mascot organizes it
                     # format 000010002000
                     variableModifications = NA,
                     # contains info on both variable and fixed modifications
                     # should be an R6 modifications class object
                     modificationTable = NA
                   ),
                   public = list(
                     #' @description
                     #' Create a new peptide object
                     #' @param sequence character vector, the amino acid
                     #'  sequence of the peptide
                     #' @param modificationTable the table from a 
                     #'  R6 'modifications' object containing the variable and
                     #'  fixed modifications present in the amino acid sequence
                     #' @param variableModifications character vector specifying
                     #'  the position of variable modifications. The length of
                     #'  this vector must be the same length as the sequence.
                     #'  Each character specifies the modification at that
                     #'  position, eg "00010", means that position 1,2,3 & 5 are
                     #'  unmodified, while position 4 has the third variable
                     #'  modification in the the modification table. Note that
                     #'  the numbering follows the original row order of the
                     #'  modification table (fixed modifications filtered out).
                     #'  Additions to a modification table should not be a
                     #'  problem, deletions or editing can cause problems
                     #'  however as the object currently cannot deal with this
                     #'  itself. If this character vector is NA, then a
                     #'  character vector of "0"'s will be created (with the
                     #'  same length as the sequence)
                     #'
                     #' @return a new 'peptide' object
                     initialize = function(sequence = "",
                                           modificationTable = NA,
                                           variableModifications= NA){
                       private$peptideString <- sequence
                       if (identical(variableModifications, NA)){
                         variableModifications <- paste(rep("0",nchar(sequence)),collapse = "")
                       }
                       if (!identical(modificationTable, NA)){
                         stopifnot(nchar(variableModifications) == nchar(sequence))
                         private$modificationTable <- modificationTable
                         private$variableModifications <- variableModifications
                       }
                       invisible(self)
                     },
                     #' @description 
                     #' For printing purposes: prints the sequence string, the
                     #'  variable modifications string and the modification
                     #'  table
                     #'
                     #' @param ... no arguments, the function takes care of
                     #'  printing
                     print = function(...){
                       cat(paste(c("Sequence               : ",self$sequence,"\n"), collapse = ""))
                       cat(paste(c("Variable modifications : ", self$modifications,"\n"), collapse = ""))
                       if (identical(self$modificationsTable, NA)){
                         cat(paste(c("Modifications Table    : ",self$modificationsTable,"\n"), collapse = ""))
                       } else {
                         cat(c("Modifications Table    :\n"))
                         print(self$modificationsTable)
                       }
                     },
                     #' @description
                     #' Retrieve part of the amino acid squence. Note: intended
                     #'  for internal use
                     #' 
                     #' @param startSeq integer vector, specifies the start of
                     #'  the part of the amino acid sequence to retrieve
                     #' @param endSeq integer vector, specifies the end of the
                     #'  part of the amino acid sequence to retrieve
                     #'
                     #' @return character vector
                     sequence.part = function(startSeq = 1L, endSeq = 1L){
                       temporary <- self$length
                       if ((startSeq %in% 1:temporary) &
                           (endSeq %in% 1:temporary) &
                           (endSeq >= startSeq)){
                         temporary <- substr(self$sequence,
                                             startSeq, endSeq)
                         return(temporary)
                       } else {
                         return(NA)
                       }
                     },
                     #' @description
                     #' Retrieve part of the variable modification string.
                     #' Note: intended for internal use
                     #' 
                     #' @param startSeq integer vector, specifies the start of
                     #'  the part of the variable modification string to
                     #'  retrieve
                     #' @param endSeq integer vector, specifies the end of the
                     #'  part of the variable modification string to retrieve
                     #'
                     #' @return character vector
                     modifications.part = function(startSeq = 1L, endSeq = 1L){
                       temporary <- self$length
                       if ((startSeq %in% 1:temporary) &
                           (endSeq %in% 1:temporary) &
                           (endSeq >= startSeq)){
                         temporary <- substr(self$modifications,
                                             startSeq, endSeq)
                         return(temporary)
                       } else {
                         return(NA)
                       }
                     },
                     #' @description 
                     #' Determines the gain & loss formulas for a part of the
                     #'  peptide (waviable modification string and modification
                     #'  table are used for this): adds up all the losses and
                     #'  gains. If the position of a variable modification in
                     #'  the variable modification string does not match the
                     #'  amino acid in the modification table, then a warning is
                     #'  produced
                     #'
                     #' @param startSeq integer vector, specifies the start of
                     #'  the part of the variable modification string to
                     #'  retrieve
                     #' @param endSeq integer vector, specifies the end of the
                     #'  part of the variable modification string to retrieve
                     #' @param Nterminal logical vector if TRUE then Nterminal
                     #'  modifications are included (if N-terminus is present in
                     #'  the part selected by startSeq and endSeq)
                     #' @param Cterminal logical vector if TRUE then Cterminal
                     #'  modifications are included (if N-terminus is present in
                     #'  the part selected by startSeq and endSeq)
                     #' 
                     #' @return a list of 2 formulas: the summed up gain
                     #'  formulas & the summed up loss formulas which are
                     #'  present in the part selected by startSeq and endSeq)
                     modifications.formula.part = function(startSeq = 1L,
                                                           endSeq = 1L,
                                                           Nterminal = TRUE,
                                                           Cterminal = TRUE){
                       temporarySequence <- self$sequence.part(startSeq,endSeq)
                       if (identical(temporarySequence,NA)){
                         return(NA)
                       } else {
                         gainFormula <- emptyFormula()
                         lossFormula <- emptyFormula()
                         temporaryMods <- self$modifications.part(startSeq,endSeq)
                         # fixed modifications
                         tempT <- self$modificationsTable %>% dplyr::filter(fixed)
                         if (nrow(tempT)>0){
                           for (counter in 1:nrow(tempT)){
                             if (tempT$position[counter] %in% c("C_term","N_term")){
                               if (tempT$position[counter] == "N_term"){
                                 if ((startSeq == 1) & Nterminal){
                                   gainFormula <- addFormulas(gainFormula, tempT$gain[[counter]])
                                   lossFormula <- addFormulas(lossFormula, tempT$loss[[counter]])
                                 } else {
                                   # no action needed
                                 }
                               } else {
                                 if ((endSeq == self$length) & Cterminal){
                                   gainFormula <- addFormulas(gainFormula, tempT$gain[[counter]])
                                   lossFormula <- addFormulas(lossFormula, tempT$loss[[counter]])
                                 } else {
                                   # no action needed
                                 }
                               }
                             } else {
                               # for every specific modification check if one of the aaa's in the position column is present
                               # note: can be more than one aa
                               positions <- unlist(strsplit(tempT$position[counter], split = ""))
                               # go over every aa in peptide string to see if it's present in positions
                               for (aa in 1:nchar(temporarySequence)){
                                 # if so, add gain and add loss from to the accumulation formulas
                                 if (substr(temporarySequence,aa,aa) %in% positions){
                                   gainFormula <- addFormulas(gainFormula, tempT$gain[[counter]])
                                   lossFormula <- addFormulas(lossFormula, tempT$loss[[counter]])
                                 }
                               }
                             }
                           }
                         }
                         # variable modifications
                         tempT <- self$modificationsTable %>% dplyr::filter(!fixed)
                         for (counter in 1:nchar(temporarySequence)){
                           currentChar <- substr(temporaryMods,counter,counter)
                           if (currentChar != "0"){
                             if (currentChar %in% 1:9){
                               varMod <- strtoi(currentChar)
                             } else {
                               # A = 10, B = 11, C = 12, etc
                               varMod <- utf8ToInt(currentChar)-55
                             }
                             gainFormula <- addFormulas(gainFormula, tempT$gain[[varMod]])
                             lossFormula <- addFormulas(lossFormula, tempT$loss[[varMod]])
                             positions <- unlist(strsplit(tempT$position[varMod], split = ""))
                             if (!(substr(temporarySequence,counter,counter) %in% positions)){
                               if (!((counter == 1) & (paste(positions, collapse = "") == "N_term"))){
                                 if (!((counter == nchar(temporarySequence)) & (paste(positions, collapse = "") == "C_term"))){
                                   warning("Variable modification positon possibly does not match the modification table!")
                                 }
                               }
                             }
                           } 
                         }
                       }
                       return(list(gain = gainFormula,loss = lossFormula))
                     },
                     #' @description 
                     #' Deterines the gain & loss formulas for the full length
                     #'  of the peptide sequence. Essentially a wrapper for
                     #'  modifications.formula.part
                     #'
                     #' @param Nterminal logical vector if TRUE then Nterminal
                     #'  modifications are included (if N-terminus is present in
                     #'  the part selected by startSeq and endSeq)
                     #' @param Cterminal logical vector if TRUE then Cterminal
                     #'  modifications are included (if N-terminus is present in
                     #'  the part selected by startSeq and endSeq)
                     #' 
                     #' @return a list of 2 formulas: the summed up gain
                     #'  formulas & the summed up loss formulas which are
                     #'  present in the part selected by startSeq and endSeq)
                     modifications.formula = function(Nterminal = TRUE, Cterminal = TRUE){
                       return(self$modifications.formula.part(startSeq = 1,
                                                           endSeq = self$length,
                                                          Nterminal = Nterminal,
                                                         Cterminal = Cterminal))
                     },
                     #' @description 
                     #' Determines the chemical formula of part of the peptide
                     #' with or without the modifications.
                     #' 
                     #' @param startSeq integer vector, specifies the start of
                     #'  the part of the peptide sequence
                     #' @param endSeq integer vector, specifies the end of the
                     #'  part of the peptide sequence
                     #' @param ignoreModifications if FALSE then modifications
                     #'  (both fixed & variable) are taken into account when
                     #'  calculating the chemical formula of the peptide. Note:
                     #'  if TRUE then the 'Nterminal' and 'Cterminal' parameters
                     #'  are ignored
                     #' @param Nterminal logical vector if TRUE then Nterminal
                     #'  modifications are included (if N-terminus is present in
                     #'  the part selected by startSeq and endSeq)
                     #' @param Cterminal logical vector if TRUE then Cterminal
                     #'  modifications are included (if N-terminus is present in
                     #'  the part selected by startSeq and endSeq)
                     #'
                     #' @return a named numeric vector, eg: c(C=6, H=12, O=6)
                     formula.part = function(startSeq = 1, endSeq = 1,
                                             ignoreModifications = FALSE,
                                            Nterminal = TRUE, Cterminal = TRUE){
                       temporary <- self$sequence.part(startSeq, endSeq)
                       #temporarymods <- self$modifications.part(startSeq, endSeq)
                       temporaryFormula <- peptideToFormula(temporary)
                       if ((!ignoreModifications) &
                           (!identical(private$modificationTable,NA))){
                    modsPresent <- self$modifications.formula.part(startSeq,
                                                                   endSeq,
                                                          Nterminal = Nterminal,
                                                          Cterminal = Cterminal)
                         temporaryFormula <- addFormulas(temporaryFormula,
                                                         modsPresent$gain)
                         temporaryFormula <- subtractFormulas(temporaryFormula,
                                                              modsPresent$loss)
                       }
                       return(temporaryFormula)
                     },
                     #' @description 
                     #' Determines the chemical formula of the full length
                     #'  peptide with or without modifications. Essentially a
                     #'  wrapper around 'formula.part'
                     #'
                     #' @param ignoreModifications if FALSE then modifications
                     #'  (both fixed & variable) are taken into account when
                     #'  calculating the chemical formula of the peptide. Note:
                     #'  if TRUE then the 'Nterminal' and 'Cterminal' parameters
                     #'  are ignored
                     #' @param Nterminal logical vector if TRUE then Nterminal
                     #'  modifications are included (if N-terminus is present in
                     #'  the part selected by startSeq and endSeq)
                     #' @param Cterminal logical vector if TRUE then Cterminal
                     #'  modifications are included (if N-terminus is present in
                     #'  the part selected by startSeq and endSeq)
                     #'
                     #' @return a named numeric vector, eg: c(C=6, H=12, O=6)
                     formula = function(ignoreModifications = FALSE,
                                        Nterminal = TRUE, Cterminal = TRUE){
                       return(self$formula.part(startSeq = 1,
                                                endSeq = self$length,
                                      ignoreModifications = ignoreModifications,
                                                Nterminal = Nterminal,
                                                Cterminal = Cterminal))
                     },
                     #' @description 
                     #' Calculate the mass of part of the peptide with or
                     #'  without modifications
                     #'
                     #' @param startSeq integer vector, specifies the start of
                     #'  the part of the peptide sequence
                     #' @param endSeq integer vector, specifies the end of the
                     #'  part of the peptide sequence
                     #' @param ignoreModifications if FALSE then modifications
                     #'  (both fixed & variable) are taken into account when
                     #'  calculating the chemical formula of the peptide. Note:
                     #'  if TRUE then the 'Nterminal' and 'Cterminal' parameters
                     #'  are ignored
                     #' @param Nterminal logical vector if TRUE then Nterminal
                     #'  modifications are included (if N-terminus is present in
                     #'  the part selected by startSeq and endSeq)
                     #' @param Cterminal logical vector if TRUE then Cterminal
                     #'  modifications are included (if N-terminus is present in
                     #'  the part selected by startSeq and endSeq)
                     #' @param elementsInfo elements masses to be used, needs to
                     #'  be of class elements, default is elementsMonoisotopic()
                     #'
                     #' @return numeric vector
                     mass.part = function(startSeq = 1, endSeq = 1,
                                          ignoreModifications = FALSE,
                                          Nterminal = TRUE, Cterminal = TRUE,
                                          elementsInfo = elementsMonoisotopic()){
                       return(self$formula.part(startSeq = startSeq,
                                                endSeq = endSeq,
                                                ignoreModifications,
                                                Nterminal = Nterminal,
                                                Cterminal = Cterminal) |>
                                formulaToMass(elementsInfo = elementsInfo))
                     },
                     #' @description 
                     #' Calculate the mass of the full length peptide with or
                     #'  without modifications
                     #'
                     #' @param ignoreModifications if FALSE then modifications
                     #'  (both fixed & variable) are taken into account when
                     #'  calculating the chemical formula of the peptide. Note:
                     #'  if TRUE then the 'Nterminal' and 'Cterminal' parameters
                     #'  are ignored
                     #' @param Nterminal logical vector if TRUE then Nterminal
                     #'  modifications are included (if N-terminus is present in
                     #'  the part selected by startSeq and endSeq)
                     #' @param Cterminal logical vector if TRUE then Cterminal
                     #'  modifications are included (if N-terminus is present in
                     #'  the part selected by startSeq and endSeq)
                     #' @param elementsInfo elements masses to be used, needs to
                     #'  be of class elements, default is elementsMonoisotopic()
                     #'
                     #' @return numeric vector
                     mass = function(ignoreModifications = FALSE,
                                     Nterminal = TRUE, Cterminal = TRUE,
                                     elementsInfo = elementsMonoisotopic()){
                       return(self$formula.part(startSeq = 1,
                                                endSeq = self$length,
                                                ignoreModifications,
                                                Nterminal = Nterminal,
                                                Cterminal = Cterminal) |>
                                formulaToMass(elementsInfo = elementsInfo))
                       
                     },
                     #' @description 
                     #' Calculate the m/z of part of the peptide (as an ion)
                     #'  with or without modifications
                     #'
                     #' @param startSeq integer vector, specifies the start of
                     #'  the part of the peptide sequence
                     #' @param endSeq integer vector, specifies the end of the
                     #'  part of the peptide sequence
                     #' @param ignoreModifications if FALSE then modifications
                     #'  (both fixed & variable) are taken into account when
                     #'  calculating the chemical formula of the peptide. Note:
                     #'  if TRUE then the 'Nterminal' and 'Cterminal' parameters
                     #'  are ignored
                     #' @param Nterminal logical vector if TRUE then Nterminal
                     #'  modifications are included (if N-terminus is present in
                     #'  the part selected by startSeq and endSeq)
                     #' @param Cterminal logical vector if TRUE then Cterminal
                     #'  modifications are included (if N-terminus is present in
                     #'  the part selected by startSeq and endSeq)
                     #' @param elementsInfo elements masses to be used, needs to
                     #'  be of class elements, default is elementsMonoisotopic()
                     #' @param adducts numeric vector, number of adducts
                     #' attached to' or 'removed from' the (originally neutral)
                     #' peptide
                     #' @param adductFormula formula (named numeric vector) of
                     #'  the adduct
                     #' @param adductCharge numeric vector indicating the actual
                     #'  charge per adduct
                     #'  
                     #' @return numeric vector
                     mz.part = function(startSeq = 1, endSeq = 1,
                                        ignoreModifications = FALSE,
                                        Nterminal = TRUE, Cterminal = TRUE,
                                        elementsInfo = elementsMonoisotopic(),
                                        adducts = 1,
                                        adductFormula = protonFormula(),
                                        adductCharge = 1){
                       return(self$mass.part(
                         startSeq = startSeq, endSeq = endSeq,
                         ignoreModifications = ignoreModifications,
                         Nterminal = Nterminal, Cterminal = Cterminal,
                         elementsInfo = elementsInfo) |>
                           massToMz(adducts = adducts,
                                    adductFormula = adductFormula,
                                    adductCharge = adductCharge,
                                    elementsInfo = elementsInfo)
                       )
                     },
                     #' @description 
                     #' Calculate the m/z of the full length peptide (as an ion)
                     #'  with or without modifications
                     #'
                     #' @param ignoreModifications if FALSE then modifications
                     #'  (both fixed & variable) are taken into account when
                     #'  calculating the chemical formula of the peptide. Note:
                     #'  if TRUE then the 'Nterminal' and 'Cterminal' parameters
                     #'  are ignored
                     #' @param Nterminal logical vector if TRUE then Nterminal
                     #'  modifications are included (if N-terminus is present in
                     #'  the part selected by startSeq and endSeq)
                     #' @param Cterminal logical vector if TRUE then Cterminal
                     #'  modifications are included (if N-terminus is present in
                     #'  the part selected by startSeq and endSeq)
                     #' @param elementsInfo elements masses to be used, needs to
                     #'  be of class elements, default is elementsMonoisotopic()
                     #' @param adducts numeric vector, number of adducts
                     #' attached to' or 'removed from' the (originally neutral)
                     #' peptide
                     #' @param adductFormula formula (named numeric vector) of
                     #'  the adduct
                     #' @param adductCharge numeric vector indicating the actual
                     #'  charge per adduct
                     #'  
                     #' @return numeric vector
                     mz = function(ignoreModifications = FALSE,
                                   Nterminal = TRUE, Cterminal = TRUE,
                                   elementsInfo = elementsMonoisotopic(),
                                   adducts = 1,
                                   adductFormula = protonFormula(),
                                   adductCharge = 1){
                       return(self$mz.part(
                         startSeq = 1, endSeq = self$length,
                         ignoreModifications = ignoreModifications,
                         Nterminal = Nterminal, Cterminal = Cterminal,
                         elementsInfo = elementsInfo,
                         adducts = adducts,
                         adductFormula = adductFormula,
                         adductCharge = adductCharge)
                       )
                     },
                     #' @description 
                     #' Calculate the m/z of part of the peptide (as a
                     #'  protonated ion) with or without modifications
                     #'  
                     #' @param startSeq integer vector, specifies the start of
                     #'  the part of the peptide sequence
                     #' @param endSeq integer vector, specifies the end of the
                     #'  part of the peptide sequence
                     #' @param ignoreModifications if FALSE then modifications
                     #'  (both fixed & variable) are taken into account when
                     #'  calculating the chemical formula of the peptide. Note:
                     #'  if TRUE then the 'Nterminal' and 'Cterminal' parameters
                     #'  are ignored
                     #' @param Nterminal logical vector if TRUE then Nterminal
                     #'  modifications are included (if N-terminus is present in
                     #'  the part selected by startSeq and endSeq)
                     #' @param Cterminal logical vector if TRUE then Cterminal
                     #'  modifications are included (if N-terminus is present in
                     #'  the part selected by startSeq and endSeq)
                     #' @param charge charge state
                     #' @param elementsInfo elements masses to be used, needs to
                     #'  be of class elements, default is elementsMonoisotopic()
                     #' 
                     #' @return numeric vector
                     mzH.part = function(startSeq = 1, endSeq = 1,
                                         ignoreModifications = FALSE,
                                         Nterminal = TRUE, Cterminal = TRUE,
                                         charge = 1,
                                         elementsInfo = elementsMonoisotopic()){
                       return(
                         self$mz.part(startSeq = startSeq,
                                      endSeq = endSeq,
                                      ignoreModifications,
                                      Nterminal = Nterminal,
                                      Cterminal = Cterminal,
                                      elementsInfo = elementsInfo,
                                      adducts = charge,
                                      adductFormula = protonFormula(),
                                      adductCharge = 1)
                       )
                     },
                     #' @description 
                     #' Calculate the m/z of part of the peptide (as a
                     #'  protonated ion) with or without modifications
                     #'  
                     #' @param charge charge state
                     #' @param ignoreModifications if FALSE then modifications
                     #'  (both fixed & variable) are taken into account when
                     #'  calculating the chemical formula of the peptide. Note:
                     #'  if TRUE then the 'Nterminal' and 'Cterminal' parameters
                     #'  are ignored
                     #' @param Nterminal logical vector if TRUE then Nterminal
                     #'  modifications are included (if N-terminus is present in
                     #'  the part selected by startSeq and endSeq)
                     #' @param Cterminal logical vector if TRUE then Cterminal
                     #'  modifications are included (if N-terminus is present in
                     #'  the part selected by startSeq and endSeq)
                     #' @param elementsInfo elements masses to be used, needs to
                     #'  be of class elements, default is elementsMonoisotopic()
                     #' 
                     #' @return numeric vector
                     mzH = function(charge = 1,
                                    ignoreModifications = FALSE,
                                    Nterminal = TRUE, Cterminal = TRUE,
                                    elementsInfo = elementsMonoisotopic()){
                       return(
                         self$mzH.part(
                           startSeq = 1, endSeq = self$length,
                           ignoreModifications,
                           Nterminal = Nterminal,
                           Cterminal = Cterminal,
                           elementsInfo = elementsInfo,
                           charge = charge)
                       )
                     },
                     #' @description generates a table of fragments which could
                     #'  arise from fragmenting part of the peptide. The
                     #'  ionseries generated are: a, a-H2O, a-NH3, b, b-H2O,
                     #'  b-NH3, b+H2O, c, x, y, y-H2O, y-NH3, z. Please note
                     #'  that the calculation is relatively 'dumb': it does NOT
                     #'  check whether a fragment is possible at all. Prime
                     #'  example is the B+H2O ion series: these fragment ions
                     #'  can only if certain conditions are met. Currently there
                     #'  is no check in this function that checks these 
                     #'  conditions/assumptions
                     #'
                     #' @param startSeq integer vector, specifies the start of
                     #'  the part of the peptide sequence
                     #' @param endSeq integer vector, specifies the end of the
                     #'  part of the peptide sequence
                     #' @param ignoreModifications if FALSE then modifications
                     #'  (both fixed & variable) are taken into account when
                     #'  calculating the chemical formula of the peptide
                     #' @param onlyIons default = TRUE, only information on the
                     #'  13 (earlier mentioned) ion series is generated. If
                     #'  FALSE then an additional 10 columns are generated with
                     #'  info on the ionseries
                     #' @param chargeState charge state of the ions in the
                     #'  generated table
                     #' @param returnFormulas default = FALSE, if TRUE then in
                     #'  stead of numerical values the table will be populated
                     #'  by the chemical formulas of the neutral fragments or
                     #'  charged fragment ions
                     #' @param formulaIncludeChargeProtons default = FALSE, if
                     #'  TRUE then protons will be included in the formulas
                     #'  (ignored when ' returnFormulas = FALSE)
                     #'
                     #' @return a data.frame with fragment information
                     fragments.part = function(startSeq = 1, endSeq = 1,
                                               ignoreModifications = FALSE,
                                               onlyIons = TRUE,
                                               chargeState = 1,
                                               returnFormulas = FALSE,
                                            formulaIncludeChargeProtons = FALSE){
                       Nterminal <- (startSeq == 1)
                       Cterminal <- (endSeq == self$length)
                       pepSeries <- list()
                       thePeptide <- private$peptideString
                       pepSeries$nTerminal <- rev(unlist(lapply((startSeq - 1):(endSeq-1),function(x){substr(thePeptide,startSeq,endSeq-x)})))
                       pepSeries$cTerminal <- rev(unlist(lapply(startSeq:endSeq,function(x){substr(thePeptide,x,endSeq)})))
                       if (!identical(private$variableModifications,NA) & !ignoreModifications){
                         thePeptide <- private$variableModifications
                         pepSeries$nVariableMods <- rev(unlist(lapply((startSeq - 1):(endSeq-1),function(x){substr(thePeptide,startSeq,endSeq-x)})))
                         pepSeries$cVariableMods <- rev(unlist(lapply(startSeq:endSeq,function(x){substr(thePeptide,x,endSeq)})))
                       } else {
                         pepSeries$nVariableMods <- rep(NA,length(pepSeries$nTerminal))
                         pepSeries$cVariableMods <- rep(NA,length(pepSeries$nTerminal))
                       }
                       pepSeries <- tibble::as_tibble(pepSeries)
                       pepSeries$nTerminalinN <- rep(Nterminal,nrow(pepSeries))
                       pepSeries$nTerminalinC <- c(rep(FALSE, nrow(pepSeries)-1),Nterminal)
                       pepSeries$cTerminalinC <- rep(Cterminal,nrow(pepSeries))
                       pepSeries$cTerminalinN <- c(rep(FALSE, nrow(pepSeries)-1),Cterminal)
                       pepSeries$nTerminalFormula <- rep(NA,nrow(pepSeries))
                       pepSeries$cTerminalFormula <- rep(NA,nrow(pepSeries))
                       for (counter in 1:nrow(pepSeries)){
                         # n-terminal
                         tempPeptide <- peptide$new(sequence = pepSeries$nTerminal[counter],
                                                    modificationTable = private$modificationTable,
                                                    variableModifications = pepSeries$nVariableMods[counter])
                         pepSeries$nTerminalFormula[counter] <- list(tempPeptide$formula(ignoreModifications = ignoreModifications,
                                                                                         Cterminal = pepSeries$cTerminalinN[counter],
                                                                                         Nterminal = pepSeries$nTerminalinN[counter]))
                         # c-terminal
                         tempPeptide <- peptide$new(sequence = pepSeries$cTerminal[counter],
                                                    modificationTable = private$modificationTable,
                                                    variableModifications = pepSeries$cVariableMods[counter])
                         pepSeries$cTerminalFormula[counter] <- list(tempPeptide$formula(ignoreModifications = ignoreModifications,
                                                                                         Cterminal = pepSeries$cTerminalinC[counter],
                                                                                         Nterminal = pepSeries$nTerminalinC[counter]))
                       }
                       if (!returnFormulas){
                         # calculate a-ions masses (not m/z's!!)
                         pepSeries$aions <- unlist(lapply(1:nrow(pepSeries),
                                                          function(x){subtractFormulas(pepSeries$nTerminalFormula[[x]],c(C=1,H=2,O=2)) %>% formulaToMass()})) %>% massToMzH(chargeState)
                         # pepSeries$aions[1] <- NA
                         # pepSeries$aions[nrow(pepSeries)] <- NA
                         # a - H2O
                         pepSeries$aionsH2O <- unlist(lapply(1:nrow(pepSeries),
                                                             function(x){subtractFormulas(pepSeries$nTerminalFormula[[x]],c(C=1,H=4,O=3)) %>% formulaToMass()})) %>% massToMzH(chargeState)
                         # a - NH3
                         pepSeries$aionsNH3 <- unlist(lapply(1:nrow(pepSeries),
                                                             function(x){subtractFormulas(pepSeries$nTerminalFormula[[x]],c(C=1,H=5,O=2,N=1)) %>% formulaToMass()})) %>% massToMzH(chargeState)                                   
                         # calculate b-ions masses (not m/z's!!)
                         pepSeries$bions <- unlist(lapply(1:nrow(pepSeries),
                                                          function(x){subtractFormulas(pepSeries$nTerminalFormula[[x]],c(H=2,O=1)) %>% formulaToMass()})) %>% massToMzH(chargeState)
                         # pepSeries$bions[1] <- NA
                         # pepSeries$bions[nrow(pepSeries)] <- NA
                         # b - H2O
                         pepSeries$bionsH2O <- unlist(lapply(1:nrow(pepSeries),
                                                             function(x){subtractFormulas(pepSeries$nTerminalFormula[[x]],c(H=4,O=2)) %>% formulaToMass()})) %>% massToMzH(chargeState)
                         # b - NH3
                         pepSeries$bionsNH3 <- unlist(lapply(1:nrow(pepSeries),
                                                             function(x){subtractFormulas(pepSeries$nTerminalFormula[[x]],c(H=5,O=1,N=1)) %>% formulaToMass()})) %>% massToMzH(chargeState)
                         # b + H2O
                         pepSeries$bionsH2Oadd <- unlist(lapply(1:nrow(pepSeries),
                                                                function(x){pepSeries$nTerminalFormula[[x]] %>% formulaToMass()})) %>% massToMzH(chargeState)
                         # calculate c-ions masses (not m/zs!!)
                         pepSeries$cions <- unlist(lapply(1:nrow(pepSeries),
                                                          function(x){
                                                            temporary <- subtractFormulas(pepSeries$nTerminalFormula[[x]],c(O=1))
                                                            return(addFormulas(temporary,c(N=1, H=1)) %>% formulaToMass())
                                                          })) %>% massToMzH(chargeState)
                         # pepSeries$cions[nrow(pepSeries)] <- NA
                         # calculate x-ions masses (not m/zs!!)
                         pepSeries$xions <- unlist(lapply(1:nrow(pepSeries),
                                                          function(x){
                                                            temporary <- subtractFormulas(pepSeries$cTerminalFormula[[x]],c(H=2))
                                                            return(addFormulas(temporary, c(C=1,O=1)) %>% formulaToMass())
                                                          })) %>% massToMzH(chargeState)
                         # pepSeries$xions[nrow(pepSeries)] <- NA
                         # calculate y-ions masses (not m/zs!!)
                         pepSeries$yions <- unlist(lapply(1:nrow(pepSeries),
                                                          function(x){pepSeries$cTerminalFormula[[x]] %>% formulaToMass()})) %>% massToMzH(chargeState)
                         # pepSeries$yions[nrow(pepSeries)] <- NA
                         # y - H2O
                         pepSeries$yionsH2O <- unlist(lapply(1:nrow(pepSeries),
                                                             function(x){subtractFormulas(pepSeries$cTerminalFormula[[x]],c(H=2,O=1)) %>% formulaToMass()})) %>% massToMzH(chargeState)
                         # y - NH3
                         pepSeries$yionsNH3 <- unlist(lapply(1:nrow(pepSeries),
                                                             function(x){subtractFormulas(pepSeries$cTerminalFormula[[x]],c(H=3,N=1)) %>% formulaToMass()})) %>% massToMzH(chargeState)
                         # calculate z-ions masses (not m/zs!!)
                         pepSeries$zions <- unlist(lapply(1:nrow(pepSeries),
                                                          function(x){subtractFormulas(pepSeries$cTerminalFormula[[x]],c(N=1,H=2)) %>% formulaToMass()})) %>% massToMzH(chargeState)
                       } else {
                         pepSeries$aions <- lapply(1:nrow(pepSeries),
                                                   function(x){subtractFormulas(pepSeries$nTerminalFormula[[x]],c(C=1,H=2,O=2))})
                         pepSeries$aionsH2O <- lapply(1:nrow(pepSeries),
                                                      function(x){subtractFormulas(pepSeries$nTerminalFormula[[x]],c(C=1,H=4,O=3))})
                         pepSeries$aionsNH3 <- lapply(1:nrow(pepSeries),
                                                      function(x){subtractFormulas(pepSeries$nTerminalFormula[[x]],c(C=1,H=5,O=2,N=1))})
                         pepSeries$bions <- lapply(1:nrow(pepSeries),
                                                   function(x){subtractFormulas(pepSeries$nTerminalFormula[[x]],c(H=2,O=1))})
                         pepSeries$bionsH2O <- lapply(1:nrow(pepSeries),
                                                      function(x){subtractFormulas(pepSeries$nTerminalFormula[[x]],c(H=4,O=2))})
                         pepSeries$bionsNH3 <- lapply(1:nrow(pepSeries),
                                                      function(x){subtractFormulas(pepSeries$nTerminalFormula[[x]],c(H=5,O=1,N=1))})
                         pepSeries$bionsH2Oadd <- lapply(1:nrow(pepSeries),
                                                         function(x){pepSeries$nTerminalFormula[[x]]})
                         pepSeries$cions <- lapply(1:nrow(pepSeries),
                                                   function(x){
                                                     temporary <- subtractFormulas(pepSeries$nTerminalFormula[[x]],c(O=1))
                                                     return(addFormulas(temporary,c(N=1, H=1)))
                                                   })
                         pepSeries$xions <- lapply(1:nrow(pepSeries),
                                                   function(x){
                                                     temporary <- subtractFormulas(pepSeries$cTerminalFormula[[x]],c(H=2))
                                                     return(addFormulas(temporary, c(C=1,O=1)))
                                                   })
                         pepSeries$yions <- lapply(1:nrow(pepSeries),
                                                   function(x){pepSeries$cTerminalFormula[[x]]})
                         pepSeries$yionsH2O <- lapply(1:nrow(pepSeries),
                                                      function(x){subtractFormulas(pepSeries$cTerminalFormula[[x]],c(H=2,O=1))})
                         pepSeries$yionsNH3 <- lapply(1:nrow(pepSeries),
                                                      function(x){subtractFormulas(pepSeries$cTerminalFormula[[x]],c(H=3,N=1))})
                         pepSeries$zions <- lapply(1:nrow(pepSeries),
                                                   function(x){subtractFormulas(pepSeries$cTerminalFormula[[x]],c(N=1,H=2))})
                         if (formulaIncludeChargeProtons){
                           for (counter in 1:ncol(pepSeries)){
                             if (grepl(colnames(pepSeries)[counter],pattern = "ions")){
                               for (rowCounter in 1:nrow(pepSeries)){
                                 pepSeries[[rowCounter, counter]][[1]] <- addFormulas(pepSeries[[rowCounter, counter]][[1]], protonFormula() * chargeState)
                               }
                             }
                           }
                         }
                       }
                       # pepSeries$zions[nrow(pepSeries)] <- NA
                       if (onlyIons) {
                         return(pepSeries[,-c(1:10)])
                       } else {
                         return(pepSeries)
                       }
                     },
                     #' @description generates a table of fragments which could
                     #'  arise from fragmenting the full sequence of the
                     #'  peptide. The ion series generated are: a, a-H2O, a-NH3,
                     #'  b, b-H2O, b-NH3, b+H2O, c, x, y, y-H2O, y-NH3, z.
                     #'  Please note that the calculation is relatively 'dumb':
                     #'  it does NOT check whether a fragment is possible at
                     #'  all. Prime example is the B+H2O ion series: these
                     #'  fragment ions can only if certain conditions are met.
                     #'  Currently there is no check in this function that
                     #'  checks these conditions/assumptions
                     #'
                     #' @param ignoreModifications if FALSE then modifications
                     #'  (both fixed & variable) are taken into account when
                     #'  calculating the chemical formula of the peptide
                     #' @param onlyIons default = TRUE, only information on the
                     #'  13 (earlier mentioned) ion series is generated. If
                     #'  FALSE then an additional 10 columns are generated with
                     #'  info on the ionseries
                     #' @param chargeState charge state of the ions in the
                     #'  generated table
                     #' @param returnFormulas default = FALSE, if TRUE then in
                     #'  stead of numerical values the table will be populated
                     #'  by the chemical formulas of the neutral fragments or
                     #'  charged fragment ions
                     #' @param formulaIncludeChargeProtons default = FALSE, if
                     #'  TRUE then protons will be included in the formulas
                     #'  (ignored when ' returnFormulas = FALSE)
                     #'
                     #' @return a data.frame with fragment information
                     fragments = function(ignoreModifications = FALSE,
                                          onlyIons = TRUE, chargeState = 1,
                                          returnFormulas = FALSE,
                                          formulaIncludeChargeProtons = FALSE){
                       self$fragments.part(startSeq = 1, endSeq = self$length, ignoreModifications = ignoreModifications,
                                           onlyIons = onlyIons, chargeState = chargeState, returnFormulas = returnFormulas, formulaIncludeChargeProtons = formulaIncludeChargeProtons)
                     },
                     #' @description generates a numeric vector containing
                     #'  'expected' immonium ions based on the amino acid
                     #'  content of part of the peptide. Please note that this
                     #'  function does NOT take into account possible
                     #'  (fixed or variable) modifications
                     #'
                     #' @param startSeq integer vector, specifies the start of
                     #'  the part of the peptide sequence
                     #' @param endSeq integer vector, specifies the end of the
                     #'  part of the peptide sequence
                     #'  
                     #' @return numeric vector
                     fragments.part.immoniumIons = function(startSeq = 1,
                                                            endSeq = 1){
                       theSequence <- self$sequence.part(startSeq,endSeq)
                       if (!is.na(theSequence)){
                         toReturn <- as.numeric()
                         theSequence <- unlist(strsplit(theSequence,split = ""))
                         for (counter in 1:length(theSequence)){
                           theIon <- switch(theSequence[counter],
                                            # in most cases:
                                            # subtractformulas(aminoAcidResidues()$getFormula("X"),c(C=1,H=0, N=0, O=1,S = 0)) %>% formulaToMass() %>% massToMzH()
                                            "D" =  88.0393,    # C3 H6 N1 O2
                                            "E" = 102.0550,    # C4 H8 N1 O2
                                            "F" = 120.0808,    # C8 H10 N1
                                            "H" =  88.0393,    # C5 H8 N3 
                                            "I" =  86.0964,    # C5 H12 N1
                                            "K" =  c(84.0808, 101.1073),    # C5 H10 N1, C5 H13 N2
                                            "L" =  86.0964,    # C5 H12 N1
                                            "M" = 104.0528,    # C4 H10 N1 S1
                                            "N" =  87.0553,    # C3 H7 N2 O1
                                            "P" = c(70.0651, 126.0550),     # C4 H8 N1, C6 H8 N1 O2 
                                            "Q" = c(84.0444, 101.0709),     # C4 H6 N1 O1, C4 H9 N2 O1
                                            "R" = c(70.0651,87.0917,100.0869,112.0869),    # C4 H8 N1, C4 H11 N2, C4 H10 N3, C5 H10 N3
                                            "S" =  60.0444,    # C2 H6 N1 O1
                                            "T" =  74.0600,    # C3 H8 N1 O1
                                            "V" =  72.0808,    # C4 H10 N1
                                            "W" = 159.0917,    # C10 H11 N2
                                            "Y" = 136.0757,    # C8 H10 N1 O1
                                            NA # default
                           )
                           if (!identical(theIon,NA)){
                             toReturn <- append(toReturn,theIon)
                           }
                         }
                         if (length(toReturn) != 0){
                           return(sort(unique(toReturn)))
                         } else {
                           return(NA)
                         }
                       } else {
                         return(NA)
                       }
                     },
                     #' @description generates a numeric vector containing
                     #'  'expected' immonium ions based on the amino acid
                     #'  content of the full sequence of the peptide. Please
                     #'  note that this function does NOT take into account
                     #'  possible (fixed or variable) modifications
                     #'
                     #' @return numeric vector
                     fragments.immoniumIons = function(){
                       self$fragments.part.immoniumIons(startSeq = 1, endSeq = self$length)
                     }
                   ),
                   active = list(
                     #' @field sequence returns the amino acid sequence as a
                     #'  character vector, can be set but is not checked against
                     #'  the length of the modifications string
                     sequence = function(value){
                       if (missing(value)){
                         return(private$peptideString)
                       } else {
                         private$peptideString <- value
                       }
                     },
                     #' @field length returns the length of the peptide (read
                     #'  only)
                     length = function(value){
                       if (missing(value)){
                         if (!identical(private$peptideString,NA)){
                           return(nchar(private$peptideString))
                         } else {
                           return(0)
                         }
                       } else {
                         # do nothing
                       }
                       
                     },
                     #' @field modifications returns the moficiations string,
                     #'  can be set but is not checked agains the length of the
                     #'  sequence string
                     modifications = function(value){
                       if (missing(value)){
                         return(private$variableModifications)
                       } else {
                         private$variableModifications <- value
                       }
                     },
                     #' @field modificationsTable returns the mofication table,
                     #'  can be modified. Note: 'variable' modifications should
                     #'  match the modifications string
                     modificationsTable = function(value){
                       if (missing(value)){
                         return(private$modificationTable)
                       } else {
                         private$modificationTable <- value
                       }
                     }
                   )
)