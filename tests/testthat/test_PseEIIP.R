context("Just testing PseEIIP functionality")

test_that("Check whether PseEIIP converts nucleotides to their electeron ion values",{
  pseEIIPfun<-as.vector(PseEIIP(seqs = "ATACG"))
  nam3mer<-nameKmer(k=3,type = "dna")
  vectTri<-vector(mode = "numeric",length = 4^3)
  names(vectTri)<-nam3mer

  EIIP<-c(0.1260,0.1335,0.1260,0.1340,0.0806)
  vectTri[c("ATA","TAC","ACG")]<-c((1/5*(0.1260+0.1335+0.1260)),(1/5*(0.1335+0.1260+0.1340)),(1/5)*(0.1260+0.1340+0.0806))
  names(vectTri)<-NULL
  expect_equal(pseEIIPfun,vectTri)
  matrix<-PseEIIP(seqs = c("ATACG","AATGG"))
  dimnames(matrix)<-NULL
  newVect<-vector(mode = "numeric",length = 4^3)
  names(newVect)<-nam3mer
  newVect[c("AAT","ATG","TGG")]<-c(1/5*(0.1260+0.1260+0.1335),1/5*(0.1260+0.1335+0.0806),1/5*(0.1335+0.0806+0.0806))
  expected_mat<-rbind(vectTri,newVect)
  dimnames(expected_mat)<-NULL
  expect_equal(matrix,expected_mat)


})
