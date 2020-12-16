context("Just testing nuc2Binary functionality")

test_that("Check whether nuc2Binary works properly",{

  diBin<-as.vector(NUC2Binary_DNA(seqs="ATAAACG",binaryType = "strBin"))
  expected<-c("0001","1000","0001","0001","0001","0010","0100")

  expect_equal(diBin,expected)



})

test_that("Check whether nuc2Binary_RNA works properly",{

  diBin<-as.vector(NUC2Binary_RNA(seqs="AUAAACG",binaryType = "strBin"))
  expected<-c("0001","1000","0001","0001","0001","0010","0100")

  expect_equal(diBin,expected)



})
