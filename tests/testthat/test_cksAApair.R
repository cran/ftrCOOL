context("Just testing CkSAAPair functionality")

test_that("Check whether CkSAAPair works with 0 sapce as the same as 2AAComposition",{
  cksaaPair<-CkSAApair(seqs="AYAAAC",rng = 0)
  expected<-kAAComposition(seqs="AYAAAC",rng = 2)
  expect_equal(cksaaPair,expected)
})

test_that("Check whether CkSAAPair works with 1 sapce as the same as 2AAComposition",{
  cksaaPair<-as.vector(CkSAApair(seqs="AYAAAC",rng = 1,normalized=FALSE))
  preVect<-c("AsA"=2,"YsA"=1,"AsC"=1)
  dipep<-nameKmer(k=2,type = "aa")

  for(i in 1:length(dipep)){
    ditemp<-unlist(strsplit(dipep[i],split = ""))
    dipep[i]<-paste(ditemp[1],ditemp[2])
  }


  featName<-vector()
  for(i in 1:1){
    featName<-c(featName,gsub(" ",strrep("s",1),dipep))
  }


  vect<-vector(length = 400,mode = "numeric")
  names(vect)<-featName
  vect[names(preVect)]<-preVect
  names(vect)<-NULL
  names(cksaaPair)<-NULL
  expect_equal(cksaaPair,vect)
})
