context("Just testing Qusian sequence order functionality")

test_that("Check whether QSOrder length vector is correct",{
  lenFun<-ncol(QSOrder(seqs = "ATSCCVYTGRILKMNSAAAFDCILG",nlag = 21))
  expectedLen<-21*2
  expect_equal(lenFun,expectedLen)
})
