context("Just testing revers Compelement functionality")

test_that("Check whether revCompm works properly",{
  revFun<-revComp(seq = "ATCGC")
  expected<-c("G","C","G","A","T")
  expect_equal(revFun,expected)
})
