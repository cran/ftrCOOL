context("Just testing SSE functionality")

test_that("Check whether length of SSEB vector is correct",{
  PredSSEdir<-system.file("testForder/",package="ftrCOOL")
  PredSSEdir<-paste0(PredSSEdir,"/test_ss2/")
  ssebFunc<-SSEB(PredSSEdir,outFormat="mat",binaryType = "numBin")
  expect_equal(dim(ssebFunc),c(2,738))
})

test_that("Check whether length of SSEC vector is correct",{
  PredSSEdir<-system.file("testForder/",package="ftrCOOL")
  PredSSEdir<-paste0(PredSSEdir,"/test_ss2/")
  ssecFunc<-SSEC(PredSSEdir)
  expect_equal(dim(ssecFunc),c(2,3))
})

test_that("Check whether length of SSES vector is correct",{
  PredSSEdir<-system.file("testForder/",package="ftrCOOL")
  PredSSEdir<-paste0(PredSSEdir,"/test_ss2/")
  sseSimpleFunc<-SSES(PredSSEdir)
  expect_equal(dim(sseSimpleFunc),c(2,246))
})

