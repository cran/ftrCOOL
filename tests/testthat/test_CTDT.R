context("Just testing CTDT functionality")

test_that("Check whether CTDT works properly",{

  ctdt<-as.vector(CTDT(seqs="AYAAACAACAFL",normalized = FALSE))
  # 222223223233
  # 232222222233
  # 131112112133

  # 222223223233
  # 131111111132
  # 212221221211
  expected9<-c(0,0,5,0,3,1,7,0,0)
  ctdt<-ctdt[1:9]
  names(expected9)<-NULL
  names(ctdt)<-NULL
  expect_equal(ctdt,expected9)

})
