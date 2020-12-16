context("Just testing maxORF functionality")

test_that("Check whether maxORF_DNA find ORF with maximum length",{

  maxorf<-as.vector(maxORF(seqs="ATGAATGCCCCGTTAACATAG"))
  expected<-"ATGAATGCCCCGTTAACATAG"
  names(expected)<-NULL
  expect_equal(maxorf,expected)
})

test_that("Check whether maxORF_DNA works properly k=1",{

  maxorf<-maxORF(seqs="ATTGAATAGCCCCGTTAACATAG")
  names(maxorf)<-NULL
  expect_equal(length(maxorf),0)
})

test_that("Check whether maxORF_RNA find ORF with maximum length",{

  maxorf<-as.vector(maxORF_RNA(seqs="AUGAAUGCCCCGUUAACAUAG"))
  expected<-"AUGAAUGCCCCGUUAACAUAG"
  names(expected)<-NULL
  expect_equal(maxorf,expected)
})
