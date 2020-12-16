context("Just testing PSTNPds functionality")

test_that("Check whether length of PSTNPds vector is correct",{

  pstFunc<-as.vector(PSTNPds(seqs="ATAAACG",pos = c("ATAAGCG","ATCCCCG"),neg = c("AAAACCG","CGCCACT")))


  expect_equal(length(pstFunc),5)

})
