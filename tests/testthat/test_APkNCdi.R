context("Just testing APkNCdi functionality")

test_that("Check whether length of APkNCdi vector is correct",{
  apknucdi<-APkNUCdi(seqs = "ACGGGCTA",selectedNucIdx = c(1,2),lambda = 3,l = 4)
  lenApknucdi<-length(apknucdi)
  expectedLen<-256+(2*3)
  expect_equal(lenApknucdi,expectedLen)
})
