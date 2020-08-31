context("Just testing Blosum62 functionality")

test_that("Check whether BLOSUM62 converts all amino acids of the sequence to a 20 dim vector ",{
  blosum_AYA<-as.vector(BLOSUM62(seqs = "AYA"))
  expected_AYA<-c(4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0,-2,-1,0
                  ,-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,
                  -1,-3,-2,-1,4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0,-2,-1,0)
  expect_equal(blosum_AYA,expected_AYA)
  blosum_AYA_ABA<-as.vector(BLOSUM62(seqs = c("AYA","ABA")))
  expect_equal(blosum_AYA_ABA,expected_AYA)

})
