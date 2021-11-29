context("Just testing KNNnucleotide functionality")

test_that("Check whether length of KNNnucleotide vector is correct",{
  knnFunc<-KNN_DNA(seqs = c("AAATCGAAT","GATCGAATC"),trainSeq = c("ACATCAACG","CTATCGAGT","AGTTCGAGT","TATTGCTAT","AAGTCGGAA","TTCTGGCAG","CCCGCTAAT","ACTTCGAAT","TCCTCGGGT","GGGTCGCAT"),percent = 20,labeltr = c(1,1,1,1,1,0,0,0,0,0))
  expect_equal(dim(knnFunc),c(2,40))
})
