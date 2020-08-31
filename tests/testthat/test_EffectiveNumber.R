context("Just testing EffectiveNumberofCodon functionality")

test_that("Check lenght of vector in Effective Number of Codon is correct",{

  LenCodonFr<-length(EffectiveNumberCodon(seqs=c("ATACGAATCATA","ATGGTCCTCATGGTGGTG")))
  expected_len<-2
  expect_equal(LenCodonFr,expected_len)

})
