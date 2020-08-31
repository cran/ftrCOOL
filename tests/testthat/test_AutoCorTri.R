context("Just testing AutoCorTriNuc functionality")

test_that("Check whether length of Moran vector is correct",{
  moranCor<-AutoCorTriNuc(seqs="AGCCCCGT",selectedNucIdx = 1:10,maxlag = 4,type = "Moran",threshold = 1)
  lenMoran<-length(moranCor)
  expected_len<-10*4
  expect_equal(lenMoran,expected_len)
})

test_that("Check whether length of Geary vector is correct",{
  gearyCor<-AutoCorTriNuc(seqs="AGCCCCGT",selectedNucIdx = 1:10,maxlag = 4,type = "Geary",threshold = 1)
  lenGeary<-length(gearyCor)
  expected_len<-10*4
  expect_equal(lenGeary,expected_len)
})

test_that("Check whether length of NormalizeMBorto vector is correct",{
  norBorCor<-AutoCorTriNuc(seqs="AGCCCCGT",selectedNucIdx = 1:10,maxlag = 4,type = "NormalizeMBorto",threshold = 1)
  lenNorBor<-length(norBorCor)
  expected_len<-10*4
  expect_equal(lenNorBor,expected_len)
})

test_that("Check whether length of AC vector is correct",{
  ACcov<-AutoCorTriNuc(seqs="AGCCCCGT",selectedNucIdx = 1:10,maxlag = 4,type = "AC",threshold = 1)
  lenAC<-length(ACcov)
  expected_len<-10*4
  expect_equal(lenAC,expected_len)
})

test_that("Check whether length of CC vector is correct",{
  CCcov<-AutoCorTriNuc(seqs="AGCCCCGT",selectedNucIdx = 1:10,maxlag = 4,type = "CC",threshold = 1)
  lenCC<-length(CCcov)
  expected_len<-(10*9*4)/2
  expect_equal(lenCC,expected_len)
})

test_that("Check whether length of ACC vector is correct",{
  ACCcov<-AutoCorTriNuc(seqs="AGCCCCGT",selectedNucIdx = 1:10,maxlag = 4,type = "ACC",threshold = 1)
  lenACC<-length(ACCcov)
  expected_len<-((10*9*4)/2)+(4*10)
  expect_equal(lenACC,expected_len)
})



