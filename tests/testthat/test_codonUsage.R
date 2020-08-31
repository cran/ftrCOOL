context("Just testing codonUsage functionality")

test_that("Check codonUsage works properly",{

  CodonUs<-CodonUsage(seqs=c("ATACGAATCATA","ATGGTCCTCATGGTGGTG"))
  expected_matrix<-matrix(0,ncol = (4^3),nrow = 2)
  colnames(expected_matrix)<-nameKmer(k=3,type = "dna")
  expected_matrix[1,c("ATA","CGA","ATC")]<-c(2,1,1)
  expected_matrix[2,c("ATG","GTC","CTC","GTG")]<-c(2,1,1,2)

  dimnames(expected_matrix)<-NULL
  dimnames(CodonUs)<-NULL
  expect_equal(CodonUs,expected_matrix)

})
