context("Just testing Dinuc2Binary functionality")

test_that("Check whether Dinuc2Binary_DNA works properly",{

  diBin<-as.vector(DiNUC2Binary_DNA(seqs="ATAAACG",binaryType = "strBin"))
  expected<-c("0011","1100","0000","0000","0001","0110")

  expect_equal(diBin,expected)



})

test_that("Check whether Dinuc2Binary_RNA works properly",{

  diBin<-as.vector(DiNUC2Binary_RNA(seqs="AUAAACG",binaryType = "strBin"))
  expected<-c("0011","1100","0000","0000","0001","0110")

  expect_equal(diBin,expected)



})

