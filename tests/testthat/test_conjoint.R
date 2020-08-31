context("Just testing conjointTriad functionality")

test_that("Check conjointTriad length is correct",{

  lenConjoint<-length(as.vector(conjointTriad(seqs="AYAAAC",normalized = FALSE)))
  expectedLen<-7*7*7
  expect_equal(lenConjoint,expectedLen)

})
