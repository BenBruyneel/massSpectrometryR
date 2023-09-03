test_that("elementsMonoisotopic works",{
  # essentially: if elementsMonoisotopic works, then elements objects work
  expect_equal(elementsMonoisotopic()$number, 17)
  expect_equal(format(elementsMonoisotopic()$getMass("Carbon"), digits = 8, nsmall = 4), "12.0000")
  expect_equal(format(elementsMonoisotopic()$getMass("O"), digits = 8, nsmall = 4), "15.994915")
  expect_equal(format(elementsMonoisotopic()$getMass("Na"), digits = 8, nsmall = 4), "22.98977")
})

test_that("elementsAverage works",{
  expect_equal(elementsAverage()$number, 13)
  expect_equal(format(elementsAverage()$getMass("Carbon"), digits = 8, nsmall = 4), "12.0107")
  expect_equal(format(elementsAverage()$getMass("O"), digits = 8, nsmall = 4), "15.9994")
  expect_equal(format(elementsAverage()$getMass("Na"), digits = 8, nsmall = 4), "22.9898")
})

test_that("sortFormula works",{
  expect_equal(sortFormula(c(C=2, H=6, O=1)), c(C=2, H=6, O=1))
  expect_equal(sortFormula(c(C=2, O=1, H=6)), c(C=2, H=6, O=1))
  expect_equal(sortFormula(c(O=1, C=2, H=6)), c(C=2, H=6, O=1))
  expect_equal(sortFormula(c(O=2, C=1, H=0)), c(C=1, H=0, O=2))
})

test_that("elementsInFormula works",{
  expect_equal(elementsInFormula(c(C=6, H=12, O=6, S=0, Na=1)),
               c("C","H","O","S","Na"))
  expect_equal(elementsInFormula(c(C=6, H=12, O=6, S=0, Na=1),
                                 removeZero = TRUE),
               c("C","H","O","Na"))
})

test_that("elementsInFormulas works",{
  expect_equal(elementsInFormulas(c(H=2,O=1), c(H=3, O= 1)),
               c("H","O"))
  expect_equal(elementsInFormulas(c(H=3,O=4, P=1), c(H=3, O= 1)),
               c("H","O","P"))
  expect_equal(elementsInFormulas(c(H=3,O=4, P=1), c(H=3, O= 1),
                                  decrease = TRUE),
               c("P","O","H"))
})

test_that("addFormulas & %f+% work",{
  expect_equal(addFormulas(c(H=2, O=1), c(H=1)),
               c(H=3, O=1))
  expect_equal(c(H=2, O=1) %f+% c(H=1),
               c(H=3, O=1))
  expect_equal(addFormulas(c(C=3, H=5, N=1, O=1, S=0), c(H=2, O=1)),
               c(C=3, H=7, N=1, O=2, S=0))
  expect_equal(c(C=3, H=5, N=1, O=1, S=0) %f+% c(H=2, O=1),
               c(C=3, H=7, N=1, O=2, S=0))
})

test_that("subtractformulas & %f-% work",{
  expect_equal(subtractFormulas(c(H=2, O=1), c(H=1)),
               c(H=1, O=1))
  expect_equal(c(H=2, O=1) %f-% c(H=1),
               c(H=1, O=1))
  expect_equal(subtractFormulas(c(C=3, H=5, N=1, O=1, S=0), c(H=2, O=1)),
               c(C=3, H=3, N=1, O=0, S=0))
  expect_equal(c(C=3, H=5, N=1, O=1, S=0) %f-% c(H=2, O=1),
               c(C=3, H=3, N=1, O=0, S=0))

})

test_that("addListFormulas works",{
  expect_equal(addListFormulas(list(c(H = 2, O = 1),
                                    c(H = 1),
                                    c(H = 2, O = 1),
                                    c(P = 1, S = 2))),
               c(H=5, O=2, P=1, S=2))
})


test_that("validFormula works",{
  expect_equal(validFormula(formula = c(C=6, H=12, O=6)), TRUE)
  expect_equal(validFormula(formula = "C6H12O6", string = T), TRUE)
  expect_equal(validFormula(formula = c(C=6, H=12, O=6, Vv = 7)), FALSE)
  expect_equal(validFormula(formula = "C6H12O6Vv7", string = T), FALSE)
})

test_that("formulaToMass works",{
  expect_equal(toString(formulaToMass(c(C=1))), "12")
  expect_equal(toString(formulaToMass(c(C=1, O=2))), "43.98982924")
  expect_equal(toString(formulaToMass(c(C = 3, H =  5, N = 1, O = 1, S = 0))),
               "71.03711379")
  expect_equal(toString(formulaToMass(c(C = 3, H =  5, N = 1, O = 1, S = 0),
                                      elementsInfo = elementsAverage())),
               "71.0777")
})

test_that("massToMz works",{
  expect_equal(format(massToMz(174.1117,
                               adducts = 1,
                               adductFormula = c(H=1),
                               adductCharge = 1),
                      digits = 8, nsmall = 4),
               "175.11898")
  expect_equal(format(massToMz(174.1117,
                               adducts = 1,
                               adductFormula = c(H=1),
                               adductCharge = 1, elementsInfo = elementsAverage()),
                      digits = 8, nsmall = 4),
               "175.1191")


})

test_that("massToMzH works",{
  expect_equal(format(massToMzH(174.1117, charge = 1),
                      digits = 8, nsmall = 4),
               "175.11898")
  expect_equal(format(massToMzH(174.1117, charge = 1,
                                elementsInfo = elementsAverage()),
                      digits = 8, nsmall = 4),
               "175.1191")
})

test_that("mzToMass works",{
  expect_equal(format(
    massToMz(mass = 174.1117, adductFormula = c(e=1), adducts = 2, adductCharge = -1) |>
      mzToMass(adductFormula = c(e=1), adducts = 2, adductCharge = -1),
    digits = 12, nsmall = 8),
    "174.11170000")
  expect_equal(format(
    massToMz(mass = 174.1117, adductFormula = c(e=1), adducts = 2, adductCharge = -1) |>
        mzToMass(adductFormula = c(e=1), adducts = 2, adductCharge = -1),
      digits = 12, nsmall = 8),
    "174.11170000")
  expect_equal(format(
    massToMz(mass = 174.1117, adductFormula = c(H=1), adducts = 2, adductCharge = 1) |>
      mzToMass(adductFormula = c(H=1), adducts = 2, adductCharge = 1),
    digits = 12, nsmall = 8),
    "174.11170000")
  expect_equal(format(
    massToMz(mass = 174.1117, adductFormula = c(H=1), adducts = 1, adductCharge = 1) |>
      mzToMass(adductFormula = c(H=1), adducts = 1, adductCharge = 1),
    digits = 12, nsmall = 8),
    "174.11170000")

})

test_that("mzHToMass works",{
  expect_equal(format(
    massToMzH(mass = 174.1117, charge = 2) |>
      mzHToMass(charge = 2),
    digits = 12, nsmall = 8),
    "174.11170000"
  )
  expect_equal(format(
    massToMzH(mass = 174.1117, charge = 1) |>
      mzHToMass(charge = 1),
    digits = 12, nsmall = 8),
    "174.11170000"
  )
})


test_that("formulaString & stringFormula work",{
  expect_equal(formulaString(c(H=3, P=1, O= 4)),
               "H3P1O4")
  expect_equal(formulaString(c(H=3, P=-1, O= 4)),
               "H3P-1O4")
  expect_equal(stringFormula("C6H12O6"), c(C=6, H=12, O=6))
  expect_equal(stringFormula("C6H12O6S-1"), c(C=6, H=12, O=6, S=-1))
})

test_that("stringToFormula works",{
  expect_equal(stringToFormula("C6H12O6"),
               c(C=6, H=12, O= 6))
  glucIsotope <- c(C=5, c13 = 1, H = 12, O18 = 1, O= 5)
  names(glucIsotope)[c(2,4)] <- c("[13]C", "[18]O")
  expect_equal(stringToFormula("C5[13]C1H12[18]O1O5"),
               glucIsotope)
})

test_that("standard formulas are correct",{
  expect_equal(protonFormula(), c(H=1))
  expect_equal(waterFormula(), c(H=2, O= 1))
  expect_equal(electronFormula(), c(e=1))
  expect_equal(emptyFormula(), c(C = 0, H =  0, N = 0, O = 0, S = 0))
})

test_that("chemicals works",{
  alcohols <- chemicals$new(
    letters = c("A","B","C"),
    shorts = c("MeOH","EthOH","PrOH"),
    names = c("Methanol","Ethanol", "Propane"),
    formulas = list(c(C=1, H=4, O=1), c(C=2, H=6, O=1), c(C=3, H=8, O=1))
  )
  expect_equal(alcohols$formulas[[1]], c(C=1, H=4, O=1))
  expect_equal(alcohols$getFormula("Methanol"), c(C=1, H=4, O=1))
  expect_equal(alcohols$getFormula("Methanl"), NA)
  expect_equal(alcohols$getFormula("MeOH"), c(C=1, H=4, O=1))
  expect_equal(alcohols$getFormula("B"), c(C=2, H=6, O=1))
  expect_equal(alcohols$letters, c("A","B","C"))
  expect_equal(alcohols$shorts, c("MeOH","EthOH","PrOH"))
  expect_equal(alcohols$names, c("Methanol","Ethanol", "Propane"))
  expect_equal(alcohols$formulas, list(c(C=1, H=4, O=1), c(C=2, H=6, O=1), c(C=3, H=8, O=1)))
})

test_that("aminoacidResidues works",{
  expect_equal(aminoAcidResidues()$number, 20)
  expect_equal(aminoAcidResidues()$getFormula("Tyrosine"), c(C=9, H=9, N=1, O=2, S =0))
  expect_equal(aminoAcidResidues()$getFormula("Tyr"), c(C=9, H=9, N=1, O=2, S =0))
  expect_equal(aminoAcidResidues()$getFormula("Y"), c(C=9, H=9, N=1, O=2, S =0))
  expect_equal(aminoAcidResidues()$getFormula("Glutamine"), c(C=5, H=8, N=2, O=2, S =0))
  expect_equal(aminoAcidResidues()$getFormula("Gln"), c(C=5, H=8, N=2, O=2, S =0))
  expect_equal(aminoAcidResidues()$getFormula("Q"), c(C=5, H=8, N=2, O=2, S =0))
})

test_that("pdToFormula works",{
  expect_equal(pdToFormula("C(6) H(12) O(6)"), c(C=6, H=12, O=6))
  expect_equal(pdToFormula("H(2) O"), c(H=2, O=1))
})

test_that("stringToFormula works",{
  expect_equal(stringToFormula("C6H12O6"), c(C=6, H=12, O=6))
  expect_equal(stringToFormula("C6H5Br"), c(C=6, H=5, Br=1))
})
