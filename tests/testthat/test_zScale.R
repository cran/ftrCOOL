context("Just testing Z_Scale functionality")

test_that("Check whether Z_Scale converts standard amino acids to zscale value",{
  zSCALE_AYA<-as.vector(zSCALE(seqs = "AYA"))
  expected_AYA<-c(0.24,-2.32,0.6,-0.14,1.3,-2.54,2.44,0.43,0.04,-1.47,0.24,-2.32,0.6,-0.14,1.3)
  expect_equal(zSCALE_AYA,expected_AYA)
  zSCALE_AYA_AAA<-zSCALE(seqs = c("AYA","AAA"))
  dimnames(zSCALE_AYA_AAA)<-NULL
  expected_AYA_AAA<-rbind(expected_AYA,c(0.24,-2.32,0.6,-0.14,1.3,0.24,-2.32,0.6,-0.14,1.3,0.24,-2.32,0.6,-0.14,1.3))
  dimnames(expected_AYA_AAA)<-NULL
  expect_equal(zSCALE_AYA_AAA,expected_AYA_AAA)

})
