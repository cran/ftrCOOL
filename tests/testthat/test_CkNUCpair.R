context("Just testing CkSNUCPair functionality")

test_that("Check whether CkSNUCPair_DNA works with 0 sapce as the same as 2NUComposition",{
  cksnucPair<-CkSNUCpair_DNA(seqs="ATCGACT",rng = 0)
  expected<-kNUComposition_DNA(seqs="ATCGACT",rng = 2)
  expect_equal(cksnucPair,expected)
})

test_that("Check whether CkSAAPair_DNA works with 1 sapce as the same as 2NUComposition",{
  cksnucPair<-as.vector(CkSNUCpair_DNA(seqs="ATCGACT",rng = 1,normalized=FALSE))
  preVect<-c("AsC"=1,"TsG"=1,"CsA"=1,"GsC"=1,"AsT"=1)
  dipep<-nameKmer(k=2,type = "dna")

  for(i in 1:length(dipep)){
    ditemp<-unlist(strsplit(dipep[i],split = ""))
    dipep[i]<-paste(ditemp[1],ditemp[2])
  }


  featName<-vector()
  for(i in 1:1){
    featName<-c(featName,gsub(" ",strrep("s",1),dipep))
  }


  vect<-vector(length = 16,mode = "numeric")
  names(vect)<-featName
  vect[names(preVect)]<-preVect
  names(vect)<-NULL
  names(cksnucPair)<-NULL
  expect_equal(cksnucPair,vect)
})

test_that("Check whether CkSNUCPair_RNA works with 0 sapce as the same as 2NUComposition",{
  cksnucPair<-CkSNUCpair_RNA(seqs="AUCGACU",rng = 0)
  expected<-kNUComposition_RNA(seqs="AUCGACU",rng = 2)
  expect_equal(cksnucPair,expected)
})
