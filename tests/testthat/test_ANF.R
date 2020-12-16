context("Just testing ANF functionality")

test_that("Check whether ANF_DNA works properly",{
  novelFun<-as.vector(ANF_DNA(seqs = "ATACG"))
  expected_AT<-c(1,1,1,1/1,0,0,1,1/2,1,1,1,2/3,0,1,0,1/4,1,0,0,1/5)
  expect_equal(novelFun,expected_AT)
  matrix_novel<-ANF_DNA(seqs = c("ATACG","AATGG"))
  dimnames(matrix_novel)<-NULL
  expected_mat<-rbind(expected_AT,c(1,1,1,1/1,1,1,1,2/2,0,0,1,1/3,1,0,0,1/4,1,0,0,2/5))
  dimnames(expected_mat)<-NULL
  expect_equal(matrix_novel,expected_mat)

  #expect_error(TriNucIndex(seqs = c("ATACGA","AAA"),selectedNucIdx = c(1,2),threshold = 1,standardized = F))

})

test_that("Check whether ANF_RNA works properly",{
  novelFun<-as.vector(ANF_RNA(seqs = "AUACG"))
  expected_AT<-c(1,1,1,1/1,0,0,1,1/2,1,1,1,2/3,0,1,0,1/4,1,0,0,1/5)
  expect_equal(novelFun,expected_AT)
  matrix_novel<-ANF_RNA(seqs = c("AUACG","AAUGG"))
  dimnames(matrix_novel)<-NULL
  expected_mat<-rbind(expected_AT,c(1,1,1,1/1,1,1,1,2/2,0,0,1,1/3,1,0,0,1/4,1,0,0,2/5))
  dimnames(expected_mat)<-NULL
  expect_equal(matrix_novel,expected_mat)

  #expect_error(TriNucIndex(seqs = c("ATACGA","AAA"),selectedNucIdx = c(1,2),threshold = 1,standardized = F))

})
