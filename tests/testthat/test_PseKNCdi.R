context("Just testing PsekNCdi functionality")

test_that("Check whether length of PsekNCdi vector is correct",{
  pseknucdi<-PSEkNCdi(seqs = "ACGGGCTA",selectedNucIdx = 1:5,lambda = 3,l = 4)
  lenpseknucdi<-length(pseknucdi)
  expectedLen<-256+(3)
  expect_equal(lenpseknucdi,expectedLen)
})
