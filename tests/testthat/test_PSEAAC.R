context("Just testing PSEAAC functionality")

test_that("Check whether length of APAAC vector is correct",{
  apaac<-PSEAAC(seqs = "AYCMLWTIL",lambda = 3,l = 2)
  lenApaac<-length(apaac)
  expectedLen<-400+(3)
  expect_equal(lenApaac,expectedLen)
})
