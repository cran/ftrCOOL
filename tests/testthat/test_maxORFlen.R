context("Just testing maxORFlength functionality")

test_that("Check whether maxORFlength_DNA find ORF with maximum length",{

  maxorf<-as.vector(maxORFlength_DNA(seqs="ATGAATGCCCCGTTAACATAG"))
  #expected<-"ATGAATGCCCCGTTAACATAG"
  expect_equal(maxorf,21)
})

test_that("Check whether maxORF works properly k=1",{

  maxorf<-maxORFlength_DNA(seqs="ATTGAATAGCCCCGTTAACATAG")
  names(maxorf)<-NULL
  expect_equal(maxorf,0)
})

test_that("Check whether maxORFlength_RNA find ORF with maximum length",{

  maxorf<-as.vector(maxORFlength_RNA(seqs="AUGAAUGCCCCGUUAACAUAG"))
  #expected<-"ATGAATGCCCCGTTAACATAG"
  expect_equal(maxorf,21)
})
