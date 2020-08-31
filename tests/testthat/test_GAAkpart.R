context("Just testing AAKpartComposition functionality")

test_that("Check AAcomposition for 2 partitions with even sequence length ",{

  comp_2part<-as.vector(GAAKpartComposition("ACYMELACYL",normalized = FALSE,k = 2,Grp = "aromatic"))
  #1 5 2 1 4 1 1 5 2 1
  excomp_p1<-rep(0,5)
  names(excomp_p1)<-nameKmer(num = 5,type = "num",k=1)
  excomp_p2<-rep(0,5)
  names(excomp_p2)<-nameKmer(num = 5,type = "num",k=1)

  Texcomp_p1<-table(unlist(strsplit("15214","")))
  Texcomp_p2<-table(unlist(strsplit("11521","")))
  excomp_p1[names(Texcomp_p1)]<-Texcomp_p1
  excomp_p2[names(Texcomp_p2)]<-Texcomp_p2

  exComp<-c(excomp_p1,excomp_p2)
  names(comp_2part)<-NULL
  names(exComp)<-NULL
  expect_equal(comp_2part,exComp)

})

