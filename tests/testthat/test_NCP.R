context("Just testing NCP functionality")

test_that("Check whether NCP_DNA works properly",{
  novelFun<-as.vector(NCP_DNA(seqs = "ATACG",binaryType = "numBin"))
  expected_AT<-c(1,1,1,0,0,1,1,1,1,0,1,0,1,0,0)
  expect_equal(novelFun,expected_AT)
  matrix_novel<-NCP_DNA(seqs = c("ATACG","AATGG"),binaryType = "numBin")
  dimnames(matrix_novel)<-NULL
  expected_mat<-rbind(expected_AT,c(1,1,1,1,1,1,0,0,1,1,0,0,1,0,0))
  dimnames(expected_mat)<-NULL
  expect_equal(matrix_novel,expected_mat)

  #expect_error(TriNucIndex(seqs = c("ATACGA","AAA"),selectedNucIdx = c(1,2),threshold = 1,standardized = F))

})

test_that("Check whether NCP_RNA works properly",{
  novelFun<-as.vector(NCP_RNA(seqs = "AUACG",binaryType = "numBin"))
  expected_AT<-c(1,1,1,0,0,1,1,1,1,0,1,0,1,0,0)
  expect_equal(novelFun,expected_AT)
  matrix_novel<-NCP_RNA(seqs = c("AUACG","AAUGG"),binaryType = "numBin")
  dimnames(matrix_novel)<-NULL
  expected_mat<-rbind(expected_AT,c(1,1,1,1,1,1,0,0,1,1,0,0,1,0,0))
  dimnames(expected_mat)<-NULL
  expect_equal(matrix_novel,expected_mat)

  #expect_error(TriNucIndex(seqs = c("ATACGA","AAA"),selectedNucIdx = c(1,2),threshold = 1,standardized = F))

})
