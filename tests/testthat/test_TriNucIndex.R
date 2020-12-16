context("Just testing TrinucIndex functionality")

test_that("Check whether TriNucIndex converts tri-nucleotides to their indexes",{
  vectTri<-as.vector(TriNUCindex_DNA(seqs = "ATACG",selectedNucIdx = c(1,2),threshold = 1))
  expected_Tri<-c(1.615,0.572,0.342,-0.07,-0.121,0.064)

  expect_equal(vectTri,expected_Tri)
  # matrix_Tri<-TriNUCindex_DNA(seqs = c("ATACG","AATGG"),selectedNucIdx = c(1,2),threshold = 1)
  # dimnames(matrix_Tri)<-NULL
  # expected_mat<-rbind(vectTri,c(0,0.35,8.7,7.7,0.7,3.05))
  # dimnames(expected_mat)<-NULL
  # expect_equal(matrix_Tri,expected_mat)

  expect_error(TriNUCindex_DNA(seqs = c("ATACGA","AAA"),selectedNucIdx = c(1,2),threshold = 1))

})
