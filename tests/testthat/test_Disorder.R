context("Just testing Disorder functionality")

test_that("Check whether length of DisorderS vector is correct",{
  PredDISdir<-system.file("testForder/",package="ftrCOOL")
  PredDISdir<-paste0(PredDISdir,"/test_Dis/")
  Disfunc<-DisorderS(PredDISdir,outFormat="mat")
  expect_equal(dim(Disfunc),c(2,246))
})

test_that("Check whether length of DisorderC vector is correct",{
  PredDISdir<-system.file("testForder/",package="ftrCOOL")
  PredDISdir<-paste0(PredDISdir,"/test_Dis/")
  Disfunc<-DisorderC(PredDISdir)
  expect_equal(dim(Disfunc),c(2,2))
})

test_that("Check whether length of DisorderB vector is correct",{
  PredDISdir<-system.file("testForder/",package="ftrCOOL")
  PredDISdir<-paste0(PredDISdir,"/test_Dis/")
  Disfunc<-DisorderB(PredDISdir,outFormat="mat",binaryType = "numBin")
  expect_equal(dim(Disfunc),c(2,492))
})
