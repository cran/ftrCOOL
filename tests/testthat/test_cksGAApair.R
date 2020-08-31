context("Just testing CkSGAApair functionality")

test_that("Check whether CkSGAAPair works with 0 sapce as the same as 2GAAComposition",{
  cksgaaPair<-CkSGAApair(seqs="AYAAAC",rng = 0,Grp = "aromatic")
  expected<-kGAAComposition(seqs="AYAAAC",rng = 2,Grp = "aromatic")
  dimnames(cksgaaPair)<-NULL
  dimnames(expected)<-NULL
  expect_equal(cksgaaPair,expected)
})

test_that("Check whether CkSGAAPair works with 1 sapce as the same as 2GAAComposition",{
  cksgaaPair<-as.vector(CkSGAApair(seqs="AYAAAC",rng = 1,normalized=FALSE,Grp = "aromatic"))
  preVect<-c("G(1s1)"=2,"G(2s1)"=1,"G(1s5)"=1)
  dipep<-nameKmer(k=2,type = "num",num = 5)

  for(i in 1:length(dipep)){
    ditemp<-unlist(strsplit(dipep[i],split = ""))
    dipep[i]<-paste(ditemp[1],ditemp[2])
  }


  featName<-vector()
  for(i in 1:1){
    featName<-c(featName,gsub(" ",strrep("s",1),dipep))
  }

  featName<-paste0("G(",featName,")")
  vect<-vector(length = 25,mode = "numeric")
  names(vect)<-featName
  vect[names(preVect)]<-preVect
  names(vect)<-NULL
  names(cksgaaPair)<-NULL
  expect_equal(cksgaaPair,vect)
})
