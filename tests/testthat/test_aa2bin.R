context("Just testing AA2Binary functionality")

test_that("Check whether AA2binary converts standard amino acids to binary value",{

  bin_AYA<-as.vector(AA2Binary("AYA",binaryType = "strBin"))
  expected_AYA<-c(paste0("1",strrep("0",19)),paste0(strrep("0",19),"1"),paste0("1",strrep("0",19)))
  expect_equal(bin_AYA,expected_AYA)
  bin_aya<-as.vector(AA2Binary("aya",binaryType = "strBin"))
  expect_equal(bin_aya,expected_AYA)
  expect_error(AA2Binary(c("AYA","AACA")))

})

test_that("testing AA2Binary with nonstandard amino acid sequences",{
  expect_message(AA2Binary(c("ABA","AYA")))
  expect_error(AA2Binary("ABA"))

})

