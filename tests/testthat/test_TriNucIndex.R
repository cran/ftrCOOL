context("Just testing TrinucIndex functionality")

test_that("Check whether TriNucIndex converts tri-nucleotides to their indexes",{
  vectTri<-as.vector(TriNucIndex(seqs = "ATACG",selectedNucIdx = c(1,2),threshold = 1,standardized = FALSE))
  expected_Tri<-c(9.7,6.25,6.4,5.05,5.2,5.3)
  expect_equal(vectTri,expected_Tri)
  matrix_Tri<-TriNucIndex(seqs = c("ATACG","AATGG"),selectedNucIdx = c(1,2),threshold = 1,standardized = FALSE)
  dimnames(matrix_Tri)<-NULL
  expected_mat<-rbind(vectTri,c(0,0.35,8.7,7.7,0.7,3.05))
  dimnames(expected_mat)<-NULL
  expect_equal(matrix_Tri,expected_mat)

  expect_error(TriNucIndex(seqs = c("ATACGA","AAA"),selectedNucIdx = c(1,2),threshold = 1,standardized = FALSE))

})
