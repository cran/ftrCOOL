context("Just testing NovelPseKNC functionality")

test_that("Check whether NovelPseKNC works properly",{
  novelFun<-as.vector(novel_PseKNC(seqs = "ATACG"))
  expected_AT<-c(1,1,0,1/1,0,0,0,1/2,1,1,0,2/3,0,1,1,1/4,1,0,1,1/5)
  expect_equal(novelFun,expected_AT)
  matrix_novel<-novel_PseKNC(seqs = c("ATACG","AATGG"))
  dimnames(matrix_novel)<-NULL
  expected_mat<-rbind(expected_AT,c(1,1,0,1/1,1,1,0,2/2,0,0,0,1/3,1,0,1,1/4,1,0,1,2/5))
  dimnames(expected_mat)<-NULL
  expect_equal(matrix_novel,expected_mat)

  #expect_error(TriNucIndex(seqs = c("ATACGA","AAA"),selectedNucIdx = c(1,2),threshold = 1,standardized = F))

})
