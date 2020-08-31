context("Just testing DDE functionality")

test_that("Check whether DDE works properly",{

  dde<-DDE(seqs="AYAAA")

  #AY=4*2   YA=2*4   AA=4*4     AA=4*4

  tm_AA=(4*4)/(61*61)
  dc_AA=2/4
  tv_AA=(tm_AA)*(1-tm_AA)/(4)
  dde_AA=(dc_AA-tm_AA)/sqrt(tv_AA)
  dde.1<-dde[1,"AA"]
  names(dde.1)<-NULL
  expect_equal(dde.1,dde_AA)

})
