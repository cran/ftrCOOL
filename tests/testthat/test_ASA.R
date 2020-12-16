context("Just testing ASA functionality")

test_that("Check whether length of ASA vector is correct",{
  PredASAdir<-system.file("testForder/",package="ftrCOOL")
  PredASAdir<-paste0(PredASAdir,"/test_ASA/")
  ASAfunc<-ASA(PredASAdir,outFormat="mat")
  expect_equal(dim(ASAfunc),c(2,246))
})
