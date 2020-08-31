context("Just testing nucKpartComposition functionality")

test_that("Check AAcomposition for 2 partitions with even sequence length ",{

  comp_2part<-as.vector(NUCKpartComposition("AGCCCCGT",normalized = FALSE,k = 2))

  excomp_p1<-rep(0,4)
  names(excomp_p1)<-c("A","C","G","T")
  excomp_p2<-rep(0,4)
  names(excomp_p2)<-c("A","C","G","T")

  Texcomp_p1<-table(unlist(strsplit("AGCC","")))
  Texcomp_p2<-table(unlist(strsplit("CCGT","")))
  excomp_p1[names(Texcomp_p1)]<-Texcomp_p1
  excomp_p2[names(Texcomp_p2)]<-Texcomp_p2

  exComp<-c(excomp_p1,excomp_p2)
  names(comp_2part)<-NULL
  names(exComp)<-NULL
  expect_equal(comp_2part,exComp)

})

test_that("functionality of AAkpartComposition for Indivisible parts ",{

  comp_2part<-as.vector(NUCKpartComposition("AGCCCCGTC",normalized = FALSE,k = 2))

  excomp_p1<-rep(0,4)
  names(excomp_p1)<-c("A","C","G","T")
  excomp_p2<-rep(0,4)
  names(excomp_p2)<-c("A","C","G","T")

  Texcomp_p1<-table(unlist(strsplit("AGCCC","")))
  Texcomp_p2<-table(unlist(strsplit("CGTC","")))
  excomp_p1[names(Texcomp_p1)]<-Texcomp_p1
  excomp_p2[names(Texcomp_p2)]<-Texcomp_p2

  exComp<-c(excomp_p1,excomp_p2)
  names(comp_2part)<-NULL
  names(exComp)<-NULL
  expect_equal(comp_2part,exComp)





  comp_3part<-as.vector(NUCKpartComposition("AGCCCCGTAGC",normalized = FALSE,k = 3))

  excomp_p1<-rep(0,4)
  names(excomp_p1)<-c("A","C","G","T")
  excomp_p2<-rep(0,4)
  names(excomp_p2)<-c("A","C","G","T")
  excomp_p3<-rep(0,4)
  names(excomp_p3)<-c("A","C","G","T")


  Texcomp_p1<-table(unlist(strsplit("AGCC","")))
  Texcomp_p2<-table(unlist(strsplit("CCGT","")))
  Texcomp_p3<-table(unlist(strsplit("AGC","")))

  excomp_p1[names(Texcomp_p1)]<-Texcomp_p1
  excomp_p2[names(Texcomp_p2)]<-Texcomp_p2
  excomp_p3[names(Texcomp_p3)]<-Texcomp_p3

  exComp<-c(excomp_p1,excomp_p2,excomp_p3)
  names(comp_3part)<-NULL
  names(exComp)<-NULL
  expect_equal(comp_3part,exComp)

})

