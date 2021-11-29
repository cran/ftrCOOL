context("Just testing ASDC functionality")

test_that("Check whether ASDC works with 0 sapce as the same as 2AAComposition",{
  cksaaPair<-as.vector(ASDC(seqs="AYAAAC"))
  names(cksaaPair)<-nameKmer(k=2,type = "aa")
  expected=vector(mode = "numeric",length = 400)
  names(expected)<-nameKmer(k=2,type = "aa")
  expected[c("AY","AA","YA","AC","YC")]<-c(1,6,3,4,1)
  sum<-15
  expected<-expected/sum
  expect_equal(cksaaPair,expected)
})

context("Just testing ASDC_DNA functionality")

test_that("Check whether ASDC_DNA works correctly or not ",{
  cksaaPair<-as.vector(ASDC_DNA(seqs="ATAAAC"))
  names(cksaaPair)<-nameKmer(k=2,type = "dna")
  expected=vector(mode = "numeric",length = 16)
  names(expected)<-nameKmer(k=2,type = "dna")
  expected[c("AT","AA","TA","AC","TC")]<-c(1,6,3,4,1)
  sum<-15
  expected<-expected/sum
  expect_equal(cksaaPair,expected)
})
