context("Just testing maxORF functionality")

test_that("Check whether maxORF find ORF with maximum length",{

  maxorf<-as.vector(maxORF(seqs="ATGAATGCCCCGTTAACATAG"))
  expected<-"ATGAATGCCCCGTTAACATAG"
  names(expected)<-NULL
  expect_equal(maxorf,expected)
})

test_that("Check whether maxORF works properly k=1",{

  maxorf<-maxORF(seqs="ATTGAATAGCCCCGTTAACATAG")
  names(maxorf)<-NULL
  expect_equal(length(maxorf),0)
})
