context("Just testing AESNN3 functionality")

test_that("Check whether AESNN3 converts standard amino acids to AESNN3 value",{
  zSCALE_AYA<-as.vector(AESNN3(seqs = "AYA"))
  expected_AYA<-c(-0.99,-0.61,0.00,0.59,0.33,-0.99,-0.99,-0.61,0.00)
  expect_equal(zSCALE_AYA,expected_AYA)
  zSCALE_AYA_AAA<-AESNN3(seqs = c("AYA","AAA"))
  dimnames(zSCALE_AYA_AAA)<-NULL
  expected_AYA_AAA<-rbind(expected_AYA,c(-0.99,-0.61,0.00,-0.99,-0.61,0.00,-0.99,-0.61,0.00))
  dimnames(expected_AYA_AAA)<-NULL
  expect_equal(zSCALE_AYA_AAA,expected_AYA_AAA)

})
