context("Just testing conjointTriadKs functionality")

test_that("Check conjointTriadKs length with 0 space with conjoint tirad",{

  ConjointKs<-conjointTriadKS(seqs="AYAAACAACAFL",rng = 0,normalized = FALSE)
  expected<-conjointTriad(seqs = "AYAAACAACAFL",normalized = FALSE)
  dimnames(ConjointKs)<-NULL
  dimnames(expected)<-NULL
  expect_equal(ConjointKs,expected)

})
