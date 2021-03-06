context("Just testing codonFraction functionality")

test_that("Check codonFraction works properly",{

  CodonFr<-CodonFraction(seqs=c("ATACGAATCATA","ATGGTCCTCATGGTGGTG"))

  codonList<-list("A"=c("GCT","GCC","GCA","GCG"),"R"=c("AGA","AGG","CGT","CGC","CGA","CGG"),
                  "N"=c("AAT","AAC"),"D"=c("GAT","GAC"),"C"=c("TGT","TGC"),"Q"=c("CAA","CAG"),
                  "E"=c("GAA","GAG"),"G"=c("GGT","GGC","GGA","GGG"),"H"=c("CAT","CAC"),
                  "I"=c("ATT","ATC","ATA"),"L"=c("TTA","TTG","CTT","CTC","CTA","CTG"),"K"=c("AAA","AAG"),
                  "M"=c("ATG"),"F"=c("TTT","TTC"),"P"=c("CCT","CCC","CCA","CCG"),
                  "S"=c("AGT","AGC","TCT","TCC","TCA","TCG"),"T"=c("ACT","ACC","ACA","ACG"),"W"=c("TGG"),
                  "Y"=c("TAT","TAC"),"V"=c("GTT","GTC","GTA","GTG"),"STOPcod"=c("TAA","TAG","TGA"))
  # I=ATA R=CGA I=ATC I=ATA
  # M=ATG V=GTC L=CTC M=ATG V=GTG v=GTG
  codonVector<-unlist(codonList)
  numCodon<-length(codonVector)
  expected_matrix<-matrix(0,ncol = numCodon,nrow = 2)
  colnames(expected_matrix)<-codonVector
  expected_matrix[1,c("ATA","CGA","ATC")]<-c(2/3,1,1/3)
  expected_matrix[2,c("ATG","GTC","CTC","GTG")]<-c(2/2,1/3,1,2/3)
  lexSort<-nameKmer(k=3,type = "dna")
  expected_matrix<-expected_matrix[,lexSort]

  dimnames(expected_matrix)<-NULL
  dimnames(CodonFr)<-NULL
  expect_equal(CodonFr,expected_matrix)

})
