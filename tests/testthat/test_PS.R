context("Just testing PS functionality")

test_that("Check whether PS2 converts standard amino acids to binary value",{

  bin_AAA<-as.vector(PS2_DNA(seqs = "AAA",binaryType = "strBin"))
  expected_AAA<-c(paste0("1",strrep("0",15)),paste0("1",strrep("0",15)))
  expect_equal(bin_AAA,expected_AAA)
  bin_aaa<-as.vector(PS2_DNA("aaa",binaryType = "strBin"))
  expect_equal(bin_aaa,expected_AAA)
  expect_error(PS2_DNA(c("AAA","AACA")))

})

test_that("testing PS2 with nonstandard amino acid sequences",{
  expect_message(PS2_DNA(c("ABA","AAA")))
  expect_error(PS2_DNA("ABA"))

})

test_that("Check whether PS3 converts standard amino acids to binary value",{

  bin_AAAA<-as.vector(PS3_DNA("AAAA",binaryType = "strBin"))
  expected_AAAA<-c(paste0("1",strrep("0",63)),paste0("1",strrep("0",63)))
  expect_equal(bin_AAAA,expected_AAAA)
  bin_aaaa<-as.vector(PS3_DNA("aaaa",binaryType = "strBin"))
  expect_equal(bin_aaaa,expected_AAAA)
  expect_error(PS3_DNA(c("AAAA","AACAA")))

})

test_that("testing PS3 with nonstandard amino acid sequences",{
  expect_message(PS3_DNA(c("ABAA","AAAA")))
  expect_error(PS3_DNA("ABAA"))

})


test_that("Check whether PS4 converts standard amino acids to binary value",{

  bin_AAAAA<-as.vector(PS4_DNA("AAAAA",binaryType = "strBin"))
  expected_AAAAA<-c(paste0("1",strrep("0",255)),paste0("1",strrep("0",255)))
  expect_equal(bin_AAAAA,expected_AAAAA)
  bin_aaaaa<-as.vector(PS4_DNA("aaaaa",binaryType = "strBin"))
  expect_equal(bin_aaaaa,expected_AAAAA)
  expect_error(PS4_DNA(c("AAAAA","AACAAA")))

})

test_that("testing PS4 with nonstandard amino acid sequences",{
  expect_message(PS4_DNA(c("ABAAA","AAAAA")))
  expect_error(PS4_DNA("ABAAA"))

})


