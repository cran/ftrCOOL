context("Just testing maxORFlength functionality")

test_that("Check whether maxORFlength find ORF with maximum length",{

  maxorf<-as.vector(maxORFlength(seqs="ATGAATGCCCCGTTAACATAG"))
  #expected<-"ATGAATGCCCCGTTAACATAG"
  expect_equal(maxorf,21)
})

test_that("Check whether maxORF works properly k=1",{

  maxorf<-maxORFlength(seqs="ATTGAATAGCCCCGTTAACATAG")
  names(maxorf)<-NULL
  expect_equal(maxorf,0)
})
