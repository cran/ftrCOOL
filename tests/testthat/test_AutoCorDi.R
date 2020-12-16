context("Just testing AutoCorDi functionality")

test_that("Check whether length of Moran_DNA vector is correct",{
  moranCor<-AutoCorDiNUC_DNA(seqs="AGCCCCGT",selectedIdx = 1:10,maxlag = 4,type = "Moran",threshold = 1)
  lenMoran<-length(moranCor)
  expected_len<-10*4
  expect_equal(lenMoran,expected_len)
})

test_that("Check whether length of Geary_DNA vector is correct",{
  gearyCor<-AutoCorDiNUC_DNA(seqs="AGCCCCGT",selectedIdx = 1:10,maxlag = 4,type = "Geary",threshold = 1)
  lenGeary<-length(gearyCor)
  expected_len<-10*4
  expect_equal(lenGeary,expected_len)
})

test_that("Check whether length of NormalizeMBorto_DNA vector is correct",{
  norBorCor<-AutoCorDiNUC_DNA(seqs="AGCCCCGT",selectedIdx = 1:10,maxlag = 4,type = "NormalizeMBorto",threshold = 1)
  lenNorBor<-length(norBorCor)
  expected_len<-10*4
  expect_equal(lenNorBor,expected_len)
})

test_that("Check whether length of AC_DNA vector is correct",{
  ACcov<-AutoCorDiNUC_DNA(seqs="AGCCCCGT",selectedIdx = 1:10,maxlag = 4,type = "AC",threshold = 1)
  lenAC<-length(ACcov)
  expected_len<-10*4
  expect_equal(lenAC,expected_len)
})

test_that("Check whether length of CC_DNA vector is correct",{
  CCcov<-AutoCorDiNUC_DNA(seqs="AGCCCCGT",selectedIdx = 1:10,maxlag = 4,type = "CC",threshold = 1)
  lenCC<-length(CCcov)
  expected_len<-(10*9*4)
  expect_equal(lenCC,expected_len)
})

test_that("Check whether length of ACC_DNA vector is correct",{
  ACCcov<-AutoCorDiNUC_DNA(seqs="AGCCCCGT",selectedIdx = 1:10,maxlag = 4,type = "ACC",threshold = 1)
  lenACC<-length(ACCcov)
  expected_len<-(10*10*4)
  expect_equal(lenACC,expected_len)
})


test_that("Check whether length of Moran_DNA vector is correct",{
  moranCor<-AutoCorDiNUC_DNA(seqs="AGCCCCGT",selectedIdx = 1:10,maxlag = 4,type = "Moran",threshold = 1)
  lenMoran<-length(moranCor)
  expected_len<-10*4
  expect_equal(lenMoran,expected_len)
})

test_that("Check whether length of Geary_DNA vector is correct",{
  gearyCor<-AutoCorDiNUC_DNA(seqs="AGCCCCGT",selectedIdx = 1:10,maxlag = 4,type = "Geary",threshold = 1)
  lenGeary<-length(gearyCor)
  expected_len<-10*4
  expect_equal(lenGeary,expected_len)
})

test_that("Check whether length of NormalizeMBorto_DNA vector is correct",{
  norBorCor<-AutoCorDiNUC_DNA(seqs="AGCCCCGT",selectedIdx = 1:10,maxlag = 4,type = "NormalizeMBorto",threshold = 1)
  lenNorBor<-length(norBorCor)
  expected_len<-10*4
  expect_equal(lenNorBor,expected_len)
})

test_that("Check whether length of AC_DNA vector is correct",{
  ACcov<-AutoCorDiNUC_DNA(seqs="AGCCCCGT",selectedIdx = 1:10,maxlag = 4,type = "AC",threshold = 1)
  lenAC<-length(ACcov)
  expected_len<-10*4
  expect_equal(lenAC,expected_len)
})

test_that("Check whether length of CC_DNA vector is correct",{
  CCcov<-AutoCorDiNUC_DNA(seqs="AGCCCCGT",selectedIdx = 1:10,maxlag = 4,type = "CC",threshold = 1)
  lenCC<-length(CCcov)
  expected_len<-(10*9*4)
  expect_equal(lenCC,expected_len)
})

test_that("Check whether length of ACC_DNA vector is correct",{
  ACCcov<-AutoCorDiNUC_DNA(seqs="AGCCCCGT",selectedIdx = 1:10,maxlag = 4,type = "ACC",threshold = 1)
  lenACC<-length(ACCcov)
  expected_len<-(10*10*4)
  expect_equal(lenACC,expected_len)
})





test_that("Check whether length of Moran_DNA vector is correct",{
  moranCor<-AutoCorDiNUC_DNA(seqs="AGCCCCGT",selectedIdx = 1:10,maxlag = 4,type = "Moran",threshold = 1)
  lenMoran<-length(moranCor)
  expected_len<-10*4
  expect_equal(lenMoran,expected_len)
})

test_that("Check whether length of Geary_DNA vector is correct",{
  gearyCor<-AutoCorDiNUC_DNA(seqs="AGCCCCGT",selectedIdx = 1:10,maxlag = 4,type = "Geary",threshold = 1)
  lenGeary<-length(gearyCor)
  expected_len<-10*4
  expect_equal(lenGeary,expected_len)
})

test_that("Check whether length of NormalizeMBorto_DNA vector is correct",{
  norBorCor<-AutoCorDiNUC_DNA(seqs="AGCCCCGT",selectedIdx = 1:10,maxlag = 4,type = "NormalizeMBorto",threshold = 1)
  lenNorBor<-length(norBorCor)
  expected_len<-10*4
  expect_equal(lenNorBor,expected_len)
})

test_that("Check whether length of AC_DNA vector is correct",{
  ACcov<-AutoCorDiNUC_DNA(seqs="AGCCCCGT",selectedIdx = 1:10,maxlag = 4,type = "AC",threshold = 1)
  lenAC<-length(ACcov)
  expected_len<-10*4
  expect_equal(lenAC,expected_len)
})

test_that("Check whether length of CC_DNA vector is correct",{
  CCcov<-AutoCorDiNUC_DNA(seqs="AGCCCCGT",selectedIdx = 1:10,maxlag = 4,type = "CC",threshold = 1)
  lenCC<-length(CCcov)
  expected_len<-(10*9*4)
  expect_equal(lenCC,expected_len)
})

test_that("Check whether length of ACC_DNA vector is correct",{
  ACCcov<-AutoCorDiNUC_DNA(seqs="AGCCCCGT",selectedIdx = 1:10,maxlag = 4,type = "ACC",threshold = 1)
  lenACC<-length(ACCcov)
  expected_len<-(10*10*4)
  expect_equal(lenACC,expected_len)
})


test_that("Check whether length of Moran_RNA vector is correct",{
  moranCor<-AutoCorDiNUC_RNA(seqs="AGCCCCGU",selectedIdx = 1:10,maxlag = 4,type = "Moran",threshold = 1)
  lenMoran<-length(moranCor)
  expected_len<-10*4
  expect_equal(lenMoran,expected_len)
})

test_that("Check whether length of Geary_RNA vector is correct",{
  gearyCor<-AutoCorDiNUC_RNA(seqs="AGCCCCGU",selectedIdx = 1:10,maxlag = 4,type = "Geary",threshold = 1)
  lenGeary<-length(gearyCor)
  expected_len<-10*4
  expect_equal(lenGeary,expected_len)
})

test_that("Check whether length of NormalizeMBorto_RNA vector is correct",{
  norBorCor<-AutoCorDiNUC_RNA(seqs="AGCCCCGU",selectedIdx = 1:10,maxlag = 4,type = "NormalizeMBorto",threshold = 1)
  lenNorBor<-length(norBorCor)
  expected_len<-10*4
  expect_equal(lenNorBor,expected_len)
})

test_that("Check whether length of AC_RNA vector is correct",{
  ACcov<-AutoCorDiNUC_RNA(seqs="AGCCCCGU",selectedIdx = 1:10,maxlag = 4,type = "AC",threshold = 1)
  lenAC<-length(ACcov)
  expected_len<-10*4
  expect_equal(lenAC,expected_len)
})

test_that("Check whether length of CC_RNA vector is correct",{
  CCcov<-AutoCorDiNUC_RNA(seqs="AGCCCCGU",selectedIdx = 1:10,maxlag = 4,type = "CC",threshold = 1)
  lenCC<-length(CCcov)
  expected_len<-(10*9*4)
  expect_equal(lenCC,expected_len)
})

test_that("Check whether length of ACC_RNA vector is correct",{
  ACCcov<-AutoCorDiNUC_RNA(seqs="AGCCCCGU",selectedIdx = 1:10,maxlag = 4,type = "ACC",threshold = 1)
  lenACC<-length(ACCcov)
  expected_len<-(10*10*4)
  expect_equal(lenACC,expected_len)
})





