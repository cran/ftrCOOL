context("Just testing EIIP functionality")

test_that("Check whether AAindex converts nucleotides to Electron-ion interaction pseudopotentials values",{
  aaindex_ATA<-as.vector(EIIP(seqs = "ATA"))
  #aaIdx<-c("A"=0.1260,"C"=0.1340,"G"=0.0806,"T"=0.1335)
  expected_ATA<-c(0.1260,0.1335,0.1260)
  expect_equal(aaindex_ATA,expected_ATA)
  aaindex_ATA_AAA<-EIIP(seqs = c("ATA","AAA"))
  dimnames(aaindex_ATA_AAA)<-NULL
  expected_ATA_AAA<-rbind(expected_ATA,c(0.1260,0.1260,0.1260))
  dimnames(expected_ATA_AAA)<-NULL
  expect_equal(aaindex_ATA_AAA,expected_ATA_AAA)

})
