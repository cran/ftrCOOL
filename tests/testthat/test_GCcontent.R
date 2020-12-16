context("Just testing G_Ccontent functionality")

test_that("Check G_Ccontent_DNA works properly",{

  gc<-G_Ccontent_DNA(seqs=c("ATACGAATCATA","ATGGTCCTCATGGTGGTG"))

  expected<-c(3/12,10/18)


  names(expected)<-NULL
  names(gc)<-NULL
  expect_equal(gc,expected)

})

test_that("Check G_Ccontent_DNA works properly",{

  gc<-G_Ccontent_RNA(seqs=c("AUACGAAUCAUA","AUGGUCCUCAUGGUGGUG"))

  expected<-c(3/12,10/18)


  names(expected)<-NULL
  names(gc)<-NULL
  expect_equal(gc,expected)

})

test_that("Check G_Ccontent_RNA works properly",{

  gc<-G_Ccontent_RNA(seqs=c("AUACGAAUCAUA","AUGGUCCUCAUGGUGGUG"))

  expected<-c(3/12,10/18)


  names(expected)<-NULL
  names(gc)<-NULL
  expect_equal(gc,expected)

})
