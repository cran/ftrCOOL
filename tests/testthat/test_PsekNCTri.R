context("Just testing PSEkNCTri functionality")

test_that("Check whether length of PSEkNCTri vector is correct",{
  pseknuctri<-PSEkNCTri(seqs = "ACGGGCTA",selectedNucIdx = c(1,2,3),lambda = 3,l = 4)
  lenpseknuctri<-length(pseknuctri)
  expectedLen<-256+(3)
  expect_equal(lenpseknuctri,expectedLen)
})
