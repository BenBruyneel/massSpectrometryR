test_that("is.Class works", {
  tester <- "Test"
  expect_equal(is.Class(tester,"character"), TRUE)
  tester <- tibble::tibble()
  expect_equal(is.Class(tester,"tbl_df"), TRUE)
  expect_equal(is.Class(tester,"tbl"), TRUE)
  expect_equal(is.Class(tester,"data.frame"), TRUE)
})

test_that("peptideFragments works",{
  expect_equal(nrow(peptideFragments()), 13)
  expect_equal(colnames(peptideFragments()), c("name","series"))
  expect_equal(peptideFragments()$name[2], "aionsH2O")
  expect_equal(peptideFragments()$series[12], "y-NH3")
})

test_that("digest works",{
  expect_equal(digest('FVNQHLCGSHLVEALYLVCGERGFFYTPKT')$peptide[2], 'GFFYTPK')
  expect_equal(digest('FVNQHLCGSHLVEALYLVCGERGFFYTPKT', missed = 2)$peptide[4],
               'FVNQHLCGSHLVEALYLVCGERGFFYTPK')
  expect_equal(nrow(digest('FVNQHLCGSHLVEALYLVCGERGFFYTPKT', missed = 2)), 6)
  expect_equal(digest('FVNQHLCGSHLVEALYLVCGERGFFYTPKT', missed = 2)$mc,
               c(0,0,0,1,1,2))
  expect_equal(nrow(digest('FVNQHLCGSHLVEALYLVCGERGFFYTPKT', enzyme = 'pepsin',
                           missed = 2)), 39)
  expect_equal(nrow(digest('FVNQHLCGSHLVEALYLVCGERGFFYTPKT',
                           enzyme = 'chymotrypsin', missed = 2)), 36)
  expect_equal(nrow(digest('FVNQHLCGSHLVEALYLVCGERGFFYTPKT',
                           enzyme = 'chymotrypsin.strict', missed = 2)), 15)
  expect_equal(nrow(digest('FVNQHLCGSHLVEALYLVCGERGFFYTPKT',
                           enzyme = 'trypsin.strict', missed = 1)), 5)
})

test_that("peptideFormula works",{
  expect_equal(peptideFormula("SAMPLER"), c(C=33, H=58, N=10, O=11, S=1))
  expect_equal(peptideFormula("SAMPLEr"), c(C=33, H=58, N=10, O=11, S=1))
  expect_equal(peptideFormula("SAMPLE"), c(C=27, H=46, N=6, O=10, S=1))
})

test_that("peptideMzH works",{
  expect_equal(format(peptideMzH("SAMPLER", charge = 1),
                      digits = 8, nsmall = 4),"803.4080")
  expect_equal(format(peptideMzH("SAMPLER", charge = 2),
               digits = 8, nsmall = 4),"402.20764")
  expect_equal(format(peptideMzH("SAMPLER", charge = 3),
               digits = 8, nsmall = 4),"268.47418")
  expect_equal(format(peptideMzH("SAMPLER",
                                 charge = 1,
                                 elementsInfo = elementsAverage()),
                      digits = 8, nsmall = 4), "803.9439")
})

test_that("peptideCount works",{
  expect_equal(peptideCount("SAMPLER", "P"), 1)
  expect_equal(peptideCount("SAMPPLER", "PLER"), 1)
  expect_equal(peptideCount("SAMPPLER", "PLER", doNotSplice = FALSE), c(2,1,1,1))
})

test_that("modifications class works (via aminoacidModifications)",{
  expect_equal(aminoAcidModifications()$number, 6)
  expect_equal(nrow(aminoAcidModifications()$fixed), 2)
  expect_equal(nrow(aminoAcidModifications()$variable), 4)
  expect_equal(addListFormulas(aminoAcidModifications()$table$gain),
               c(C=4, H=8, N=1, O=8))
  expect_equal(addListFormulas(aminoAcidModifications()$table$loss),
               c(C=0, H=4, N=1, O=0, S=0))
})

test_that("peptide class works",{
  tester <- peptide$new(sequence = "SAMPLER",
                        variableModifications = "0010000",
                        modificationTable = aminoAcidModifications()$table[-2,])
  expect_equal(tester$sequence, "SAMPLER")
  expect_equal(tester$modifications, "0010000")
  expect_equal(tester$formula(),
               c(C=33, H=58, N=10, O=12, S=1))
  expect_equal(tester$formula(ignoreModifications = TRUE),
               c(C=33, H=58, N=10, O=11, S=1))
  expect_equal(format(tester$mass(), digits= 8, nsmall = 4),
               "818.39564")
  expect_equal(format(tester$mz(adductFormula = c(Na=1)),
                      digits = 8, nsmall = 4),
               "841.38486")
  expect_equal(format(tester$mzH(), digits = 8, nsmall = 4),
               "819.40291")
  expect_equal(format(tester$mzH(charge = 2,
                                 ignoreModifications = TRUE),
                      digits = 8, nsmall = 4),
               "402.20764")
  expect_equal(format(sum(tester$fragments()$bions),
                      digits = 8, nsmall = 4),
               "2919.3244")
  expect_equal(format(sum(tester$fragments(ignoreModifications = TRUE)$bions),
                      digits = 8, nsmall = 4),
               "2839.3498")
  expect_equal(format(sum(tester$fragments()$aions),
                      digits = 8, nsmall = 4),
               "2723.3600")
  expect_equal(format(sum(tester$fragments()$aionsH2O),
                      digits = 8, nsmall = 4),
               "2597.2860")
  expect_equal(format(sum(tester$fragments()$aionsNH3),
                      digits = 8, nsmall = 4),
               "2604.1741")
  expect_equal(format(sum(tester$fragments()$bions),
                      digits = 8, nsmall = 4),
               "2919.3244")
  expect_equal(format(sum(tester$fragments()$bionsH2O),
                      digits = 8, nsmall = 4),
               "2793.2504")
  expect_equal(format(sum(tester$fragments()$bionsNH3),
                      digits = 8, nsmall = 4),
               "2800.1385")
  expect_equal(format(sum(tester$fragments()$bionsH2Oadd),
                      digits = 8, nsmall = 4),
               "3045.3983")
  expect_equal(format(sum(tester$fragments()$cions),
                      digits = 8, nsmall = 4),
               "3038.5102")
  expect_equal(format(sum(tester$fragments()$xions),
                      digits = 8, nsmall = 4),
               "3805.7869")
  expect_equal(format(sum(tester$fragments()$yions),
                      digits = 8, nsmall = 4),
               "3623.9321")
  expect_equal(format(sum(tester$fragments()$yionsH2O),
                      digits = 8, nsmall = 4),
               "3497.8581")
  expect_equal(format(sum(tester$fragments()$yionsNH3),
                      digits = 8, nsmall = 4),
               "3504.7462")
  expect_equal(format(sum(tester$fragments()$zions),
                      digits = 8, nsmall = 4),
               "3511.8010")
  expect_equal(addListFormulas(tester$fragments(returnFormulas = TRUE)$yions),
               c(C=146, H=265, N=49,O=51,S=3))
  expect_equal(tester$modifications.formula(),
               list(gain = c(C=0, H=0, N= 0, O=1, S=0),
                    loss = c(C=0, H=0, N= 0, O=0, S=0)))
  expect_equal(format(sum(tester$fragments.immoniumIons()),
                      digits = 8, nsmall = 4),
               "847.6342")
  expect_equal(nrow(tester$modificationsTable), 5)
  tester <- peptide$new(sequence = "CAMPLER",
                        variableModifications = "0010000",
                        modificationTable = aminoAcidModifications()$table[-2,])
  expect_equal(tester$formula(),
               c(C=35, H=61, N=11, O=12, S=2))
  expect_equal(tester$formula(ignoreModifications = TRUE),
               c(C=33, H=58, N=10, O=10, S=2))
  expect_equal(format(tester$mass(), digits= 8, nsmall = 4),
               "891.39426")
  expect_equal(addListFormulas(tester$fragments(returnFormulas = TRUE)$yions),
               c(C=148, H=268, N=50,O=51,S=4))
})

test_that("aminoAcidResidues works",{
  expect_equal(aminoAcidResidues()$translatePeptide(
    sequence = "SAMPLER", from1to3 = TRUE,
    joinCharacter = "-"),
    "Ser-Ala-Met-Pro-Leu-Glu-Arg")
  expect_equal(aminoAcidResidues()$translatePeptide(
    sequence = "Ser-Ala-Met-Pro-Leu-Glu-Arg", from1to3 = FALSE,
    splitCharacter = "-"),
    "SAMPLER")
})
