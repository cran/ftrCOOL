context("Just testing Dinuc2Binary functionality")

test_that("Check whether Dinuc2Binary works properly",{

  diBin<-as.vector(NUC2Binary(seqs="ATAAACG",binaryType = "strBin"))
  expected<-c("0001","1000","0001","0001","0001","0010","0100")

  expect_equal(diBin,expected)



})
