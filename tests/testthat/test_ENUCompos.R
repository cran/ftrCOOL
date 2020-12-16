context("Just testing ENUCompositon functionality")

test_that("Check whether ENUCompositon_DNA without overLAp works properly",{
  enhancedNucCompos<-as.vector(ENUComposition_DNA(seqs = "ATAATCGCC",winSize = 3,overLap = FALSE))
  nuc<-c("A","C","G","T")
  vect1<-vector(length = 4,mode = "numeric")
  names(vect1)<-nuc
  vect1["A"]=2
  vect1["T"]=1

  vect2<-vector(length = 4,mode = "numeric")
  names(vect2)<-nuc
  vect2["A"]=1
  vect2["T"]=1
  vect2["C"]=1

  vect3<-vector(length = 4,mode = "numeric")
  names(vect3)<-nuc
  vect3["G"]=1
  vect3["C"]=2

  expected<-c(vect1,vect2,vect3)
  names(expected)<-NULL
  expect_equal(enhancedNucCompos,expected)


})

test_that("Check whether ENUCompositon_DNA with overLAp works properly",{
  enhancedNucCompos<-as.vector(ENUComposition_DNA(seqs = "ATAATCGCC",winSize = 7,overLap = TRUE))
  nuc<-c("A","C","G","T")
  vect1<-vector(length = 4,mode = "numeric")
  names(vect1)<-nuc
  vect1["A"]=3
  vect1["T"]=2
  vect1["C"]=1
  vect1["G"]=1

  vect2<-vector(length = 4,mode = "numeric")
  names(vect2)<-nuc
  vect2["A"]=2
  vect2["T"]=2
  vect2["C"]=2
  vect2["G"]=1


  vect3<-vector(length = 4,mode = "numeric")
  names(vect3)<-nuc
  vect3["A"]=2
  vect3["T"]=1
  vect3["C"]=3
  vect3["G"]=1


  expected<-c(vect1,vect2,vect3)
  names(expected)<-NULL
  expect_equal(enhancedNucCompos,expected)


})

test_that("Check whether ENUCompositon_DNA without overLAp works properly",{
  enhancedNucCompos<-as.vector(ENUComposition_DNA(seqs = "ATAATCGC",winSize = 3,overLap = FALSE))
  nuc<-c("A","C","G","T")
  vect1<-vector(length = 4,mode = "numeric")
  names(vect1)<-nuc
  vect1["A"]=2
  vect1["T"]=1

  vect2<-vector(length = 4,mode = "numeric")
  names(vect2)<-nuc
  vect2["A"]=1
  vect2["T"]=1
  vect2["C"]=1

  vect3<-vector(length = 4,mode = "numeric")
  names(vect3)<-nuc
  vect3["G"]=1
  vect3["C"]=1

  expected<-c(vect1,vect2,vect3)
  names(expected)<-NULL
  expect_equal(enhancedNucCompos,expected)


})

test_that("Check whether ENUCompositon_RNA without overLAp works properly",{
  enhancedNucCompos<-as.vector(ENUComposition_RNA(seqs = "AUAAUCGCC",winSize = 3,overLap = FALSE))
  nuc<-c("A","C","G","U")
  vect1<-vector(length = 4,mode = "numeric")
  names(vect1)<-nuc
  vect1["A"]=2
  vect1["U"]=1

  vect2<-vector(length = 4,mode = "numeric")
  names(vect2)<-nuc
  vect2["A"]=1
  vect2["U"]=1
  vect2["C"]=1

  vect3<-vector(length = 4,mode = "numeric")
  names(vect3)<-nuc
  vect3["G"]=1
  vect3["C"]=2

  expected<-c(vect1,vect2,vect3)
  names(expected)<-NULL
  expect_equal(enhancedNucCompos,expected)


})

