context("Just testing PsekNCdi functionality")

test_that("Check whether length of PsekNCdi_DNA vector is correct",{
  pseknucdi<-PSEkNUCdi_DNA(seqs = "ACGGGCTA",selectedIdx = 1:5,lambda = 3,l = 4)
  lenpseknucdi<-length(pseknucdi)
  expectedLen<-256+(3)
  expect_equal(lenpseknucdi,expectedLen)
})

test_that("Check whether length of PsekNCdi_RNA vector is correct",{
  pseknucdi<-PSEkNUCdi_RNA(seqs = "ACGGGCUA",selectedIdx = 1:5,lambda = 3,l = 4)
  lenpseknucdi<-length(pseknucdi)
  expectedLen<-256+(3)
  expect_equal(lenpseknucdi,expectedLen)
})
