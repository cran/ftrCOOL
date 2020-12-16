context("Just testing Torsion functionality")

test_that("Check whether length of Torsion vector is correct",{
  PredTordir<-system.file("testForder/",package="ftrCOOL")
  PredTordir<-paste0(PredTordir,"/test_Torsion/")
  torFunc<-TorsionAngle(PredTordir,outFormat="mat")
  expect_equal(dim(torFunc),c(2,492))
})
