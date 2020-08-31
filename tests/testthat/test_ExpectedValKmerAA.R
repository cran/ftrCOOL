context("Just testing Expected Value kmer amino acids functionality")

test_that("Check whether ExpectedValKmerAA works properly",{
  valKmer<-as.vector(ExpectedValueKmerAA(seqs="ACTCMTC",k = 2,normalized = FALSE))
  #AT TC CA AT
  vect=vector(mode = "numeric",length = 400)
  names(vect)=nameKmer(k=2,type = "aa")
  vect["AC"]=1/(1*3)
  vect["CT"]=1/(3*2)
  vect["TC"]=2/(2*3)
  vect["CM"]=1/(3*1)
  vect["MT"]=1/(1*2)
  names(vect)<-NULL
  expect_equal(valKmer,vect)
})

