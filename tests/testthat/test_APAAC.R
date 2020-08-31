context("Just testing APAAC functionality")

test_that("Check whether length of APAAC vector is correct",{
  apaac<-APAAC(seqs = "AYCMLWTIL",lambda = 3,l = 2)
  lenApaac<-length(apaac)
  expectedLen<-400+(2*3)
  expect_equal(lenApaac,expectedLen)
})
