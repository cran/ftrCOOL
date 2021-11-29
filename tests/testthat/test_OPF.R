context("Just testing OPF functionality")

test_that("Check whether OPF_10bit works properly",{

  opf<-as.vector(OPF_10bit(seqs="AYAAA"))
  A=c(0,0,0,0,1,0,1,0,1,0)
  Y=c(1,0,0,1,1,0,0,0,0,0)
  expected<-c(A,Y,A,A,A)
  names(expected)<-NULL
  names(opf)<-NULL
  expect_equal(opf,expected)

})


test_that("Check whether OPF_7bit works properly",{

  opf<-as.vector(OPF_7bit_T1(seqs="AYAAA"))
  A=c(1,0,1,0,1,0,1)
  Y=c(1,0,0,1,0,0,0)
  expected<-c(A,Y,A,A,A)
  names(expected)<-NULL
  names(opf)<-NULL
  expect_equal(opf,expected)

})
