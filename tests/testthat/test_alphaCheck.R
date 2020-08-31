context("Just testing alphabet_Check functionality")

test_that("Check whether nonStandard amino acid alphabet",{
  expect_error(alphabetCheck(sequences = "AAABA",alphabet = "aa"))
  expect_message(alphabetCheck(sequences = c("AAABA","AAAACA"),alphabet = "aa"))

})

test_that("Check whether DNA alphabet is correct",{
  expect_error(alphabetCheck(sequences = "UGAC",alphabet = "dna"))
  expect_message(alphabetCheck(sequences = c("UGAC","ATCGA"),alphabet = "dna"))

})

test_that("Check whether RNA alphabet is correct",{
  expect_error(alphabetCheck(sequences = "TGAC",alphabet = "rna"))
  expect_message(alphabetCheck(sequences = c("TGAC","AUCGA"),alphabet = "rna"))

})

