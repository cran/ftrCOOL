context("Just testing localPoSpkNUCF functionality")

test_that("Check whether localPoSpkNUCF works properly(Vector)",{
  localspeFreq<-as.vector(LocalPoSpKNUCF_DNA(seqs="AAGAGCC",k=2))
  expected<-c(1/2,1/3,1/4,2/5,1/6,1/7)
  expect_equal(localspeFreq,expected)
})

test_that("Check whether localPoSpkaaF works properly(Matrix)",{
  localspeFreq<-LocalPoSpKNUCF_DNA(seqs=c("AAGAGCC","ACCCACC"),k=2)
  expected<-rbind(c(1/2,1/3,1/4,2/5,1/6,1/7),c(1/2,1/3,2/4,1/5,2/6,3/7))
  dimnames(localspeFreq)<-NULL
  dimnames(expected)<-NULL
  expect_equal(localspeFreq,expected)
})

test_that("Check localPoSpkaaF for sequences with different length",{
  expect_error(LocalPoSpKNUCF_DNA(seqs=c("AAGAGCC","ACC"),k=2))
})

test_that("Check whether localPoSpkNUCF works properly(Vector)",{
  localspeFreq<-as.vector(LocalPoSpKNUCF_RNA(seqs="AAGAGCC",k=2))
  expected<-c(1/2,1/3,1/4,2/5,1/6,1/7)
  expect_equal(localspeFreq,expected)
})
