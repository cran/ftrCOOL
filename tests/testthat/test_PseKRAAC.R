context("Just testing PseKRAAC functionality")

test_that("Check whether PseKRAAC_T1 with lambda type works properly",{

  psekrracLam<-as.vector(PseKRAAC_T1(seqs="AYAAAC",type = "lambda",Grp = 4,k = 2,GapOrLambdaValue = 2))
  names<-nameKmer(type = "num",k = 2,num = 4)
  featurVect<-vector(mode = "numeric",length = length(names))
  names(featurVect)<-names
  featurVect["33"]=1
  featurVect["13"]=1
  featurVect["31"]=1
  names(featurVect)<-NULL
  expect_equal(psekrracLam,featurVect)

})

test_that("Check whether PseKRAAC_T1 with gap type works properly",{

  psekrracGap<-as.vector(PseKRAAC_T1(seqs="AYAAAC",type = "gap",Grp = 4,k = 2,GapOrLambdaValue = 2))
  names<-nameKmer(type = "num",k = 2,num = 4)
  featurVectEx<-vector(mode = "numeric",length = length(names))
  names(featurVectEx)<-names
  featurVectEx["31"]=2
  names(featurVectEx)<-NULL
  expect_equal(psekrracGap,featurVectEx)

})


