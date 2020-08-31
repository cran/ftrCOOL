context("Just testing CTDC functionality")

test_that("Check whether CTDC works properly",{

  ctdc<-as.vector(CTDC(seqs="AYAAACAACAFL",normalized = FALSE))
  # 222223223233
  # 232222222233
  # 131112112133
  expected9<-c(0,8,4,0,9,3,7,2,3)
  ctdc<-ctdc[1:9]
  names(expected9)<-NULL
  names(ctdc)<-NULL
  expect_equal(ctdc,expected9)

})
