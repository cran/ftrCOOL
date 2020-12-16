context("Just testing APkNCdi functionality")

test_that("Check whether length of APkNCdi_DNA vector is correct",{
  apknucdi<-APkNUCdi_DNA(seqs = "ACGGGCTA",selectedIdx = c(1,2),lambda = 3,l = 4)
  lenApknucdi<-length(apknucdi)
  expectedLen<-256+(2*3)
  expect_equal(lenApknucdi,expectedLen)
})

test_that("Check whether length of APkNCdi_RNA vector is correct",{
  apknucdi<-APkNUCdi_RNA(seqs = "ACGGGCUA",selectedIdx = c(1,2),lambda = 3,l = 4)
  lenApknucdi<-length(apknucdi)
  expectedLen<-256+(2*3)
  expect_equal(lenApknucdi,expectedLen)
})
