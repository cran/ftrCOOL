context("Just testing Dinuc2Binary functionality")

test_that("Check whether Dinuc2Binary works properly",{

  diBin<-as.vector(Dinuc2Binary(seqs="ATAAACG",binaryType = "strBin"))
  expected<-c("0011","1100","0000","0000","0001","0110")

  expect_equal(diBin,expected)



})
