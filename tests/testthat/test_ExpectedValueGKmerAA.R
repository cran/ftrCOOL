context("Just testing Expected Value Grouped kmer amino acids functionality")

test_that("Check whether ExpectedValGKmerAA works properly",{
  valKmer<-as.vector(ExpectedValueGKmerAA(seqs="ACTCMTC",k = 2,normalized = FALSE,Grp = "aromatic"))
  #AT TC CA AT
    #  Group =list(Grp1=c("G","A","V","L","M","I")
  #,Grp2=c("F","Y","W")(2,2,1)=5,Grp3=c("K","R","H")(2,6,2)=10
  #,Grp4=c("D","E")(2,2)=4,Grp5=c("S","T","C","P","N","Q"))(6,4,2,4,2,2)=20

  vect=vector(mode = "numeric",length = 25)
  names(vect)=nameKmer(k=2,type = "num",num = 5)
  vect["15"]=2/(2*5)
  vect["55"]=3/(5*5)
  vect["51"]=1/(5*2)

  names(vect)<-NULL
  expect_equal(valKmer,vect)
})

