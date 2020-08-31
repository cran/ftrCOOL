context("Just testing Expected Value Grouped amino acids functionality")

test_that("Check whether ExpectedValGAA works properly",{
  valKmer<-as.vector(ExpectedValueGAA(seqs="ACTCMTC",k = 2,normalized = FALSE,Grp = "aromatic"))
  #AT TC CA AT
  numCodAA <- list("A"=4,"C"=2,"D"=2,"E"=2,"F"=2,"G"=4,"H"=2,"I"=3,"K"=2,"L"=6,"M"=1,"N"=2,"P"=4,"Q"=2,"R"=6,"S"=6,"T"=4,"V"=4,"W"=1,"Y"=2)
  #  Group =list(Grp1=c("G","A","V","L","M","I")(4,4,4,6,1,3)=22
  #,Grp2=c("F","Y","W")(2,2,1)=5,Grp3=c("K","R","H")(2,6,2)=10
  #,Grp4=c("D","E")(2,2)=4,Grp5=c("S","T","C","P","N","Q"))(6,4,2,4,2,2)=20

  vect=vector(mode = "numeric",length = 25)
  names(vect)=nameKmer(k=2,type = "num",num = 5)
  vect["15"]=2/(22*20)
  vect["55"]=3/(20*20)
  vect["51"]=1/(20*22)

  names(vect)<-NULL
  expect_equal(valKmer,vect)
})

