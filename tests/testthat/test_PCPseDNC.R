context("Just testing PCPseDNC functionality")

test_that("Check whether length of PCPseDNC vector is correct",{
  pseknucdi<-PCPseDNC(seqs = "ACGGGCTA",lambda = 3,l = 4)
  lenpseknucdi<-length(pseknucdi)
  expectedLen<-256+(3)
  expect_equal(lenpseknucdi,expectedLen)
})
