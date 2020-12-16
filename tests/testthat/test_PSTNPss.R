context("Just testing PSTNPss functionality")

test_that("Check whether length of PSTNPss_DNA vector is correct",{

  pstFunc<-as.vector(PSTNPss_DNA(seqs="ATAAACG",pos = c("ATAAGCG","ATCCCCG"),neg = c("AAAACCG","CGCCACT")))


  expect_equal(length(pstFunc),5)

})

test_that("Check whether length of PSTNPss_RNA vector is correct",{

  pstFunc<-as.vector(PSTNPss_RNA(seqs="AUAAACG",pos = c("AUAAGCG","AUCCCCG"),neg = c("AAAACCG","CGCCACU")))


  expect_equal(length(pstFunc),5)

})
