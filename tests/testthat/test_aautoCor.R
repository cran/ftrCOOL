context("Just testing AAutoCor functionality")

test_that("Check whether length of Moran vector is correct",{
  moranCor<-AAutoCor(seqs="AYAAAC",maxlag = 3,type = "Moran",threshold = 1)
  lenMoran<-length(moranCor)
  expected_len<-24
  expect_equal(lenMoran,expected_len)
})

test_that("Check whether length of Geary vector is correct",{
  gearyCor<-AAutoCor(seqs="AYAAAC",maxlag = 3,type = "Geary",threshold = 1)
  lenGeary<-length(gearyCor)
  expected_len<-24
  expect_equal(lenGeary,expected_len)
})

test_that("Check whether length of NormalizeMBorto vector is correct",{
  norBorCor<-AAutoCor(seqs="AYAAAC",maxlag = 3,type = "NormalizeMBorto",threshold = 1)
  lenNorBor<-length(norBorCor)
  expected_len<-24
  expect_equal(lenNorBor,expected_len)
})

test_that("Check whether length of AC vector is correct",{
  ACcov<-AAutoCor(seqs="AYAAAC",maxlag = 3,type = "AC",threshold = 1)
  lenAC<-length(ACcov)
  expected_len<-24
  expect_equal(lenAC,expected_len)
})

test_that("Check whether length of CC vector is correct",{
  CCcov<-AAutoCor(seqs="AYAAAC",maxlag = 3,type = "CC",threshold = 1)
  lenCC<-length(CCcov)
  expected_len<-8*7*3
  expect_equal(lenCC,expected_len)
})

test_that("Check whether length of ACC vector is correct",{
  ACCcov<-AAutoCor(seqs="AYAAAC",maxlag = 3,type = "ACC",threshold = 1)
  lenACC<-length(ACCcov)
  expected_len<-8*8*3
  expect_equal(lenACC,expected_len)
})



