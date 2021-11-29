context("Just testing binary_bit functionality")

test_that("Check whether binary6bit converts standard amino acids to binary value",{

  bin_AYA<-as.vector(binary_6bit("AYA",binaryType  = "numBin"))
  expected_AYA<-c(0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,1,0,0)
  expect_equal(bin_AYA,expected_AYA)
  bin_aya<-as.vector(binary_6bit("aya",binaryType = "numBin"))
  expect_equal(bin_aya,expected_AYA)
  expect_error(binary_6bit(c("AYA","AACA")))

})

test_that("testing binary6bit with nonstandard amino acid sequences",{
  expect_message(binary_6bit(c("ABA","AYA")))
  expect_error(binary_6bit("ABA"))

})

test_that("Check whether binary3bit converts standard amino acids to binary value",{

  bin_AYA<-as.vector(binary_3bit_T1("AYA",binaryType  = "numBin"))
  expected_AYA<-c(0,1,0,0,1,0,0,1,0)
  expect_equal(bin_AYA,expected_AYA)
  bin_aya<-as.vector(binary_3bit_T1("aya",binaryType = "numBin"))
  expect_equal(bin_aya,expected_AYA)
  expect_error(binary_3bit_T1(c("AYA","AACA")))

})

test_that("Check whether binary5bit converts standard amino acids to binary value",{

  bin_AYA<-as.vector(binary_5bit_T1("AYA",binaryType  = "numBin"))
  expected_AYA<-c(1,0,0,0,0,0,1,0,0,0,1,0,0,0,0)
  expect_equal(bin_AYA,expected_AYA)
  bin_aya<-as.vector(binary_5bit_T1("aya",binaryType = "numBin"))
  expect_equal(bin_aya,expected_AYA)
  expect_error(binary_5bit_T1(c("AYA","AACA")))

})


