context("Just testing EAACompositon functionality")

test_that("Check whether EAACompositon without overLAp works properly",{
  enhancedAACompos<-as.vector(EAAComposition(seqs = "AYAATCGLC",winSize = 3,overLap = FALSE))
  aa<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
  vect1<-vector(length = 20,mode = "numeric")
  names(vect1)<-aa
  vect1["A"]=2
  vect1["Y"]=1

  vect2<-vector(length = 20,mode = "numeric")
  names(vect2)<-aa
  vect2["A"]=1
  vect2["T"]=1
  vect2["C"]=1

  vect3<-vector(length = 20,mode = "numeric")
  names(vect3)<-aa
  vect3["G"]=1
  vect3["L"]=1
  vect3["C"]=1

  expected<-c(vect1,vect2,vect3)
  names(expected)<-NULL
  expect_equal(enhancedAACompos,expected)


})

test_that("Check whether EAACompositon with overLAp works properly",{
  enhancedAACompos<-as.vector(EAAComposition(seqs = "AYAATCGLC",winSize = 7,overLap = TRUE))
  aa<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
  vect1<-vector(length = 20,mode = "numeric")
  names(vect1)<-aa
  vect1["A"]=3
  vect1["Y"]=1
  vect1["T"]=1
  vect1["C"]=1
  vect1["G"]=1

  vect2<-vector(length = 20,mode = "numeric")
  names(vect2)<-aa
  vect2["A"]=2
  vect2["Y"]=1
  vect2["T"]=1
  vect2["C"]=1
  vect2["G"]=1
  vect2["L"]=1

  vect3<-vector(length = 20,mode = "numeric")
  names(vect3)<-aa
  vect3["A"]=2
  vect3["T"]=1
  vect3["C"]=2
  vect3["G"]=1
  vect3["L"]=1

  expected<-c(vect1,vect2,vect3)
  names(expected)<-NULL
  expect_equal(enhancedAACompos,expected)


})

test_that("Check whether EAACompositon without overLAp works properly",{
  enhancedAACompos<-as.vector(EAAComposition(seqs = "AYAATCGL",winSize = 3,overLap = FALSE))
  aa<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
  vect1<-vector(length = 20,mode = "numeric")
  names(vect1)<-aa
  vect1["A"]=2
  vect1["Y"]=1

  vect2<-vector(length = 20,mode = "numeric")
  names(vect2)<-aa
  vect2["A"]=1
  vect2["T"]=1
  vect2["C"]=1

  vect3<-vector(length = 20,mode = "numeric")
  names(vect3)<-aa
  vect3["G"]=1
  vect3["L"]=1

  expected<-c(vect1,vect2,vect3)
  names(expected)<-NULL
  expect_equal(enhancedAACompos,expected)


})
