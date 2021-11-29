context("Just testing Mismatch functionality")

test_that("Check whether Mismatch works properly",{
  mismatch<-Mismatch_DNA(seqs="ATCGACT",k = 2,m=1)

  expectedAA<-3
  expect_equal(mismatch[1,"AA"],expectedAA)
})
