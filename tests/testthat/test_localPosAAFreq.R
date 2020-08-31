context("Just testing localPoSpkaaF functionality")

test_that("Check whether localPoSpkaaF works properly(Vector)",{
  localspeFreq<-as.vector(LocalPoSpKaaF(seqs="AAMAMCI",k=2))
  expected<-c(1/2,1/3,1/4,2/5,1/6,1/7)
  expect_equal(localspeFreq,expected)
})

test_that("Check whether localPoSpkaaF works properly(Matrix)",{
  localspeFreq<-LocalPoSpKaaF(seqs=c("AAMAMCI","AMMMAMM"),k=2)
  expected<-rbind(c(1/2,1/3,1/4,2/5,1/6,1/7),c(1/2,1/3,2/4,1/5,2/6,3/7))
  dimnames(localspeFreq)<-NULL
  dimnames(expected)<-NULL
  expect_equal(localspeFreq,expected)
})

test_that("Check localPoSpkaaF for sequences with different length",{
  expect_error(LocalPoSpKaaF(seqs=c("AAMAMCI","AMM"),k=2))
})
