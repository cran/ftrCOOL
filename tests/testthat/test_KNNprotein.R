context("Just testing KNNprotein functionality")

test_that("Check whether length of KNNprotein vector is correct",{
  knnFunc<-KNNProtein(seqs = c("MARQHPSAL","GARQHPSAL"),trainSeq = c("MARQHPSAA","MARQHHSAS","AARQHPSAA","AARRHPSSE","AARLHTSSE","AYRIHPTVE","AYRIHPTVE","QYIIHPTVE","MARQHPSAS","MARQHPSAL"),percent = 20,labeltr = c(1,1,1,1,1,0,0,0,0,0))
  expect_equal(dim(knnFunc),c(2,40))
  # checkmat=matrix(0,ncol = 40,nrow = 2)
  # a1=(0:9)*2+1
  # a2=(10:19)*2+1
  # checkmat[,a1]=1
  # checkmat[,a2]=2
  # colnames(knnFunc)<-NULL
  # row.names(knnFunc)<-NULL
  # colnames(checkmat)<-NULL
  # row.names(checkmat)<-NULL
  # knnFunc<-as.matrix(knnFunc)
  # expect_equal(knnFunc,checkmat)

})
