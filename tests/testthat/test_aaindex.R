context("Just testing AAindex functionality")

test_that("Check whether AAindex converts standard amino acids to physicochemical index",{
  aaindex_AYA<-as.vector(AAindex(seqs = "AYA",selectedAAidx = c(1,2),threshold = 1,standardized = FALSE))
  expected_AYA<-c(4.35,0.61,4.6,1.88,4.35,0.61)
  expect_equal(aaindex_AYA,expected_AYA)
  aaindex_AYA_AAA<-AAindex(seqs = c("AYA","AAA"),selectedAAidx = c(1,2),threshold = 1,standardized = FALSE)
  dimnames(aaindex_AYA_AAA)<-NULL
  expected_AYA_AAA<-rbind(expected_AYA,c(4.35,0.61,4.35,0.61,4.35,0.61))
  dimnames(expected_AYA_AAA)<-NULL
  expect_equal(aaindex_AYA_AAA,expected_AYA_AAA)

})
