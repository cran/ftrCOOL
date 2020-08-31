context("Just testing DiNucIndex functionality")

test_that("Check whether DiNucIndex converts di-nucleotides to their indexes",{
  vectDi<-as.vector(DiNucIndex(seqs = "ATACG",selectedNucIdx = c(1,2),threshold = 1,standardized = FALSE))
  expected_Di<-c(33.81,-10.6,33.28,-11.2,31.12,-11.8,32.91,-13.1)
  expect_equal(vectDi,expected_Di)
  matrix_Di<-DiNucIndex(seqs = c("ATACG","AATGG"),selectedNucIdx = c(1,2),threshold = 1,standardized = FALSE)
  dimnames(matrix_Di)<-NULL
  expected_mat<-rbind(vectDi,c(38.9,-12,33.81,-10.6,41.41,-12.3,34.96,-9.5))
  dimnames(expected_mat)<-NULL
  expect_equal(matrix_Di,expected_mat)

  expect_error(DiNucIndex(seqs = c("ATACGA","AAA"),selectedNucIdx = c(1,2),threshold = 1,standardized = FALSE))

})
