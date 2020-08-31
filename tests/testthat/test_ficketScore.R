context("Just testing ficketScore functionality")

test_that("Check whether ficketScore is a number",{
  ficket<-length(fickettScore(seqs="ATCGACT"))
  expected<-1
  expect_equal(ficket,expected)
})

