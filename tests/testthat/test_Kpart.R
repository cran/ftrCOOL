context("Just testing AAKpartComposition functionality")

test_that("Check AAcomposition for 2 partitions with even sequence length ",{

  comp_2part<-as.vector(AAKpartComposition("ACYMELACYL",normalized = FALSE,k = 2))

  excomp_p1<-rep(0,20)
  names(excomp_p1)<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y" )
  excomp_p2<-rep(0,20)
  names(excomp_p2)<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y" )

  Texcomp_p1<-table(unlist(strsplit("ACYME","")))
  Texcomp_p2<-table(unlist(strsplit("LACYL","")))
  excomp_p1[names(Texcomp_p1)]<-Texcomp_p1
  excomp_p2[names(Texcomp_p2)]<-Texcomp_p2

  exComp<-c(excomp_p1,excomp_p2)
  names(comp_2part)<-NULL
  names(exComp)<-NULL
  expect_equal(comp_2part,exComp)

})

test_that("functionality of AAkpartComposition for Indivisible parts ",{

  comp_2part<-as.vector(AAKpartComposition("ACYMELACY",normalized = F,k = 2))

  excomp_p1<-rep(0,20)
  names(excomp_p1)<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y" )
  excomp_p2<-rep(0,20)
  names(excomp_p2)<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y" )

  Texcomp_p1<-table(unlist(strsplit("ACYME","")))
  Texcomp_p2<-table(unlist(strsplit("LACY","")))
  excomp_p1[names(Texcomp_p1)]<-Texcomp_p1
  excomp_p2[names(Texcomp_p2)]<-Texcomp_p2

  exComp<-c(excomp_p1,excomp_p2)
  names(comp_2part)<-NULL
  names(exComp)<-NULL
  expect_equal(comp_2part,exComp)





  comp_3part<-as.vector(AAKpartComposition("ACYMELACYME",normalized = F,k = 3))

  excomp_p1<-rep(0,20)
  names(excomp_p1)<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y" )
  excomp_p2<-rep(0,20)
  names(excomp_p2)<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y" )
  excomp_p3<-rep(0,20)
  names(excomp_p3)<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y" )


  Texcomp_p1<-table(unlist(strsplit("ACYM","")))
  Texcomp_p2<-table(unlist(strsplit("ELAC","")))
  Texcomp_p3<-table(unlist(strsplit("YME","")))

  excomp_p1[names(Texcomp_p1)]<-Texcomp_p1
  excomp_p2[names(Texcomp_p2)]<-Texcomp_p2
  excomp_p3[names(Texcomp_p3)]<-Texcomp_p3

  exComp<-c(excomp_p1,excomp_p2,excomp_p3)
  names(comp_3part)<-NULL
  names(exComp)<-NULL
  expect_equal(comp_3part,exComp)

})

