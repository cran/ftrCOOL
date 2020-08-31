context("Just testing sequence order coupling number functionality")

test_that("Check whether SOCNumber length vector is correct",{
  lenFun<-ncol(SOCNumber(seqs = "ATSCCVYTGRILKMNSAAAFDCILG",nlag = 21))
  expectedLen<-21*2
  expect_equal(lenFun,expectedLen)
})
