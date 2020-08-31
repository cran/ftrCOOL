context("Just testing CTD functionality")

test_that("Check whether CTD works properly",{

  ctd<-as.vector(CTD(seqs="AYAAACAACAFL",normalized=FALSE))
  # 222223223233
  # 232222222233
  # 131112112133
  expectedc<-as.vector(CTDC(seqs="AYAAACAACAFL",normalized = FALSE))
  expectedt<-as.vector(CTDT(seqs="AYAAACAACAFL",normalized = FALSE))
  expectedd<-as.vector(CTDD(seqs="AYAAACAACAFL"))

  expectedctd<-c(expectedc,expectedt,expectedd)


  names(ctd)<-NULL
  names(expectedctd)<-NULL
  expect_equal(ctd,expectedctd)

})
