context("Just testing CTDD functionality")

test_that("Check whether CTDD works properly",{

  ctdd<-as.vector(CTDD(seqs="AYAAACAACAFL"))
  # 222223223233
  # 232222222233
  # 131112112133
  expected15<-c(0,0,0,0,0,1/12,2/12,4/12,7/12,10/12,6/12,6/12,9/12,11/12,12/12)
  ctdd<-ctdd[1:15]
  names(expected15)<-NULL
  names(ctdd)<-NULL
  expect_equal(ctdd,expected15)

})
