context("Just testing APkNCTri functionality")

test_that("Check whether length of APkNCTri vector is correct",{
  apknuctri<-APkNUCTri(seqs = "ACGGGCTA",selectedNucIdx = c(1,2,3),lambda = 3,l = 4)
  lenApknuctri<-length(apknuctri)
  expectedLen<-256+(3*3)
  expect_equal(lenApknuctri,expectedLen)
})
