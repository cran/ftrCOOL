context("Just testing Expected Value amino acids functionality")

test_that("Check whether ExpectedValAA works properly",{
  valKmer<-as.vector(ExpectedValueAA(seqs="ACTCMTC",k = 2,normalized = FALSE))
  #AT TC CA AT
  numCodAA <- list("A"=4,"C"=2,"D"=2,"E"=2,"F"=2,"G"=4,"H"=2,"I"=3,"K"=2,"L"=6,"M"=1,"N"=2,"P"=4,"Q"=2,"R"=6,"S"=6,"T"=4,"V"=4,"W"=1,"Y"=2)
  vect=vector(mode = "numeric",length = 400)
  names(vect)=nameKmer(k=2,type = "aa")
  vect["AC"]=1/(4*2)
  vect["CT"]=1/(2*4)
  vect["TC"]=2/(4*2)
  vect["CM"]=1/(2*1)
  vect["MT"]=1/(1*4)
  names(vect)<-NULL
  expect_equal(valKmer,vect)
})

