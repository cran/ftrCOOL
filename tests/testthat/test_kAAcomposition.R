context("Just testing kAAComposition functionality")

test_that("Check whether AAComposition works properly k=1",{

  kAAcompos<-as.vector(kAAComposition(seqs="AYAAAC",rng = 1,normalized = FALSE))
  expected<-vector(mode = "numeric",length = 20)
  names(expected)<-nameKmer(k=1,type = "aa")
  expected[c("A","Y","C")]=c(4,1,1)
  names(expected)<-NULL
  expect_equal(kAAcompos,expected)
})

test_that("Check whether Dipeptide Composition works properly",{
  kAAcompos<-as.vector(kAAComposition(seqs="AYAAAC",rng = 2,normalized = FALSE))
  expected<-vector(mode = "numeric",length = 400)
  names(expected)<-nameKmer(k=2,type = "aa")
  expected[c("AY","YA","AA","AC")]=c(1,1,2,1)
  names(expected)<-NULL
  expect_equal(kAAcompos,expected)
})
