context("Just testing MMI functionality")

test_that("Check whether length of MMI vector is correct",{
  apaac<-MMI_DNA(seqs = "ACCGACGT")
  lenApaac<-length(apaac)
  expectedLen<-30
  expect_equal(lenApaac,expectedLen)
})
