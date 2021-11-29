context("Just testing DistancePair functionality")

test_that("Check whether DistancePair works with 0 sapce as the same as 2GAAComposition",{
  cksgaaPair<-DistancePair(seqs="AYAAAC",rng = 0,Grp = "cp14")
  expected<-kGAAComposition(seqs="AYAAAC",rng = 2,Grp = list(Grp1=c("M","I","V")
                                                                    ,Grp2=c("L"),Grp3=c("F"),Grp4=c("W","Y"),Grp5=c("G")
                                                                    ,Grp6=c("P"),Grp7=c("C"),Grp8=c("A"),Grp9=c("S"),Grp10=c("T"),Grp11=c("N"),Grp12=c("H","R","K","Q"),Grp13=c("E"),Grp14=c("D")))

  expectedAC<-kGAAComposition(seqs="AYAAAC",rng = 1,Grp = list(Grp1=c("M","I","V")
                                                             ,Grp2=c("L"),Grp3=c("F"),Grp4=c("W","Y"),Grp5=c("G")
                                                             ,Grp6=c("P"),Grp7=c("C"),Grp8=c("A"),Grp9=c("S"),Grp10=c("T"),Grp11=c("N"),Grp12=c("H","R","K","Q"),Grp13=c("E"),Grp14=c("D")))

  expected<-cbind(expectedAC,expected)
  dimnames(cksgaaPair)<-NULL
  dimnames(expected)<-NULL
  expect_equal(cksgaaPair,expected)
})
