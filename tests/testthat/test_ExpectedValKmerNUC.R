context("Just testing Expected Value kmer nucleotide functionality")

test_that("Check whether ExpectedValKmerNUC_DNA works properly",{
  valKmer<-as.vector(ExpectedValKmerNUC_DNA(seqs="ATCAT",k = 2,normalized = FALSE))
  #AT TC CA AT
  vect=vector(mode = "numeric",length = 16)
  names(vect)=nameKmer(k=2,type = "dna")
  vect["AT"]=2/(2*2)
  vect["TC"]=1/(2*1)
  vect["CA"]=1/(1*2)
  names(vect)<-NULL
  expect_equal(valKmer,vect)
})

test_that("Check whether ExpectedValKmerNUC_RNA works properly",{
  valKmer<-as.vector(ExpectedValKmerNUC_RNA(seqs="AUCAU",k = 2,normalized = FALSE))
  #AT TC CA AT
  vect=vector(mode = "numeric",length = 16)
  names(vect)=nameKmer(k=2,type = "rna")
  vect["AU"]=2/(2*2)
  vect["UC"]=1/(2*1)
  vect["CA"]=1/(1*2)
  names(vect)<-NULL
  expect_equal(valKmer,vect)
})


