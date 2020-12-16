context("Just testing PSSM functionality")

test_that("Check whether length of PSSM vector is correct",{
  PredPSSMdir<-system.file("testForder/",package="ftrCOOL")
  PredPSSMdir<-paste0(PredPSSMdir,"/test_PSSM/")
  pssmFunc<-PSSM(PredPSSMdir,outFormat="mat")
  expect_equal(dim(pssmFunc),c(2,4920))
})
