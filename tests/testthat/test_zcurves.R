context("Just testing Zcurve functionality")

test_that("Check whether length of Zcurve9bit vector is correct",{
  apaac<-Zcurve9bit_DNA(seqs = "ACCGACGT")
  lenApaac<-length(apaac)
  expectedLen<-9
  expect_equal(lenApaac,expectedLen)
})

test_that("Check whether length of Zcurve12bit vector is correct",{
  apaac<-Zcurve12bit_DNA(seqs = "ACCGACGT")
  lenApaac<-length(apaac)
  expectedLen<-12
  expect_equal(lenApaac,expectedLen)
})

test_that("Check whether length of Zcurve36bit vector is correct",{
  apaac<-Zcurve36bit_DNA(seqs = "ACCGACGT")
  lenApaac<-length(apaac)
  expectedLen<-36
  expect_equal(lenApaac,expectedLen)
})

test_that("Check whether length of Zcurve48bit vector is correct",{
  apaac<-Zcurve48bit_DNA(seqs = "ACCGACGT")
  lenApaac<-length(apaac)
  expectedLen<-48
  expect_equal(lenApaac,expectedLen)
})

test_that("Check whether length of Zcurve144bit vector is correct",{
  apaac<-Zcurve144bit_DNA(seqs = "ACCGACGT")
  lenApaac<-length(apaac)
  expectedLen<-144
  expect_equal(lenApaac,expectedLen)
})


