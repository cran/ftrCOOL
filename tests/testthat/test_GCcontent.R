context("Just testing codonFraction functionality")

test_that("Check codonFraction works properly",{

  gc<-GCcontent(seqs=c("ATACGAATCATA","ATGGTCCTCATGGTGGTG"))

  expected<-c(3/12,10/18)


  names(expected)<-NULL
  names(gc)<-NULL
  expect_equal(gc,expected)

})
