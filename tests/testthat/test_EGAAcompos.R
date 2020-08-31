context("Just testing EAACompositon functionality")

test_that("Check whether EAACompositon without overLAp works properly",{
  enhancedGAACompos<-as.vector(EGAAComposition(seqs = "AYAATCGLC",winSize = 3,overLap = FALSE,Grp = "aromatic"))
  aa<-c("1","2","3","4","5")
  #A=1 Y=2 T=5 C=5 G=1 L=1
  vect1<-vector(length = 5,mode = "numeric")
  names(vect1)<-aa
  vect1["1"]=2
  vect1["2"]=1

  vect2<-vector(length = 5,mode = "numeric")
  names(vect2)<-aa
  vect2["1"]=1
  vect2["5"]=2


  vect3<-vector(length = 5,mode = "numeric")
  names(vect3)<-aa
  vect3["1"]=2
  vect3["5"]=1

  expected<-c(vect1,vect2,vect3)
  names(expected)<-NULL
  expect_equal(enhancedGAACompos,expected)


})

test_that("Check whether EAACompositon with overLap works properly",{
  enhancedGAACompos<-as.vector(EGAAComposition(seqs = "AYAATCGLC",winSize = 7,overLap = TRUE,Grp="aromatic"))
  #A=1 Y=2 T=5 C=5 G=1 L=1
  aa<-c("1","2","3","4","5")
  vect1<-vector(length = 5,mode = "numeric")
  names(vect1)<-aa
  vect1["1"]=4
  vect1["2"]=1
  vect1["5"]=2

  vect2<-vector(length = 5,mode = "numeric")
  names(vect2)<-aa
  vect2["1"]=4
  vect2["2"]=1
  vect2["5"]=2

  vect3<-vector(length = 5,mode = "numeric")
  names(vect3)<-aa
  vect3["1"]=4
  vect3["5"]=3

  expected<-c(vect1,vect2,vect3)
  names(expected)<-NULL
  expect_equal(enhancedGAACompos,expected)


})

test_that("Check whether EAACompositon without overLAp works properly",{
  enhancedGAACompos<-as.vector(EGAAComposition(seqs = "AYAATCGL",winSize = 3,overLap = FALSE,Grp="aromatic"))
  #A=1 Y=2 T=5 C=5 G=1 L=1
  aa<-c("1","2","3","4","5")
  vect1<-vector(length = 5,mode = "numeric")
  names(vect1)<-aa
  vect1["1"]=2
  vect1["2"]=1

  vect2<-vector(length = 5,mode = "numeric")
  names(vect2)<-aa
  vect2["1"]=1
  vect2["5"]=2

  vect3<-vector(length = 5,mode = "numeric")
  names(vect3)<-aa
  vect3["1"]=2

  expected<-c(vect1,vect2,vect3)
  names(expected)<-NULL
  expect_equal(enhancedGAACompos,expected)


})
