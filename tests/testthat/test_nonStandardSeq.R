context("Just testing nonStandardSeq functionality")

test_that("Check whether nonStandardSeq works properly",{

  seq1="ABCAAMN"
  seq2="AERTYLIQ"
  seq3="HCVLIYTME"
  seq4="SAALJUZX"

  seqs1<-c("seq1"=seq1,"seq2"=seq2,"seq3"=seq3)
  seqs2<-c("s1"=seq1,"s2"=seq2,"s4"=seq4)
  nonStan1<-nonStandardSeq(seqs1)
  expect_equal(nonStan1,seqs1[1])

  nonStan2<-nonStandardSeq(seqs2)
  expect_equal(nonStan2,c(seqs2[1],seqs2[3]))
})
