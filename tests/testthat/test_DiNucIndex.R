context("Just testing DiNucIndex functionality")

test_that("Check whether DiNucIndex_ converts di-nucleotides to their indexes",{
  vectDi<-as.vector(DiNUCindex_DNA(seqs = "ATACG",selectedIdx = c(1,2),threshold = 1))
  #expected_Di<-c(33.81,-10.6,33.28,-11.2,31.12,-11.8,32.91,-13.1)
  expected_Di<-c(0.567,-1.05,1.603,0.418,-0.918,-0.831,-0.579,2.229)
  expect_equal(vectDi,expected_Di)
  # matrix_Di<-DiNucIndex(seqs = c("ATACG","AATGG"),selectedIdx = c(1,2),threshold = 1)
  # dimnames(matrix_Di)<-NULL
  # expected_mat<-rbind(vectDi,c(38.9,-12,33.81,-10.6,41.41,-12.3,34.96,-9.5))
  # dimnames(expected_mat)<-NULL
  # expect_equal(matrix_Di,expected_mat)

  expect_error(DiNUCindex_DNA(seqs = c("ATACGA","AAA"),selectedIdx = c(1,2),threshold = 1))

})

test_that("Check whether DiNucIndex converts di-nucleotides to their indexes",{
  vectDi<-as.vector(DiNUCindex_RNA(seqs = "AUACG",selectedIdx = c(1,2),threshold = 1))
  expected_Di<-c(-1.36,1,-1.45,1,-1.43,1,-1.89,0)
  expect_equal(vectDi,expected_Di)
  # matrix_Di<-DiNucIndex(seqs = c("ATACG","AATGG"),selectedNucIdx = c(1,2),threshold = 1)
  # dimnames(matrix_Di)<-NULL
  # expected_mat<-rbind(vectDi,c(38.9,-12,33.81,-10.6,41.41,-12.3,34.96,-9.5))
  # dimnames(expected_mat)<-NULL
  #expect_equal(matrix_Di,expected_mat)

  expect_error(DiNUCindex_RNA(seqs = c("AUACGA","AAA"),selectedNucIdx = c(1,2),threshold = 1))

})
