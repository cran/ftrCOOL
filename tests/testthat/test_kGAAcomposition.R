context("Just testing kGAAComposition functionality")

test_that("Check whether kGAAComposition works properly k=1",{

  kGAAcompos<-as.vector(kGAAComposition(seqs="AYAAACV",rng = 1,normalized = FALSE,Grp = "aromatic"))
  #1211151
  expected<-vector(mode = "numeric",length = 5)
  names(expected)<-nameKmer(k=1,type = "num",num = 5)
  expected[c("1","2","5")]=c(5,1,1)
  names(expected)<-NULL
  expect_equal(kGAAcompos,expected)
})

test_that("Check whether Dipeptide Composition works properly",{
  kGAAcompos<-as.vector(kGAAComposition(seqs="AYAAACV",rng = 2,normalized = FALSE,Grp = "aromatic"))
  expected<-vector(mode = "numeric",length = 25)
  names(expected)<-nameKmer(k=2,type = "num",num = 5)
  expected[c("12","21","11","15","51")]=c(1,1,2,1,1)
  names(expected)<-NULL
  expect_equal(kGAAcompos,expected)
})
