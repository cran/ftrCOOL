context("Just testing needleman functionality")

test_that("Check Effective needleman works properly",{
  seq1="HELLO_WORLD"
  seq2="HELLOWORLD"
  need=needleman(seq1,seq2)
  expect_equal(need,9)

})
