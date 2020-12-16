#' Pseudo Electron-Ion Interaction Pseudopotentials of Trinucleotide (PseEIIP)
#'
#' This function calculates the pseudo electron-ion interaction for each sequence.
#' It creates a feature vector for each sequence.
#' The vector contains a value for each for each tri-nucleotide.
#' The value is computed by multiplying the aggregate value of electron-ion interaction of each nucleotide
#'
#'
#'
#' @references Chen, Zhen, et al. "iLearn: an integrated platform and meta-learner for feature engineering, machine-learning analysis and modeling of DNA, RNA and protein sequence data." Briefings in bioinformatics 21.3 (2020): 1047-1057.
#'
#'
#' @param seqs is a FASTA file containing nucleotide sequences. The sequences start
#' with '>'. Also, seqs could be a string vector. Each element of the vector is a nucleotide sequence.
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#'
#' @return This function returns a feature matrix which the number of rows is equal to the number of sequences and the
#' number of columns is 4^3=64.
#'
#'
#' @export
#'
#' @examples
#'
#' LNCSeqsADR<-system.file("extdata/",package="ftrCOOL")
#' LNC50Nuc<-as.vector(read.csv(paste0(LNCSeqsADR,"/LNC50Nuc.csv"))[,2])
#' mat<-PseEIIP(seqs = LNC50Nuc)


PseEIIP <- function(seqs,label=c())
{

  if(length(seqs)==1&&file.exists(seqs)){
    seqs<-fa.read(seqs,alphabet="dna")
    seqs_Lab<-alphabetCheck(seqs,alphabet = "dna",label)

    seqs<-seqs_Lab[[1]]
    label<-seqs_Lab[[2]]
  }
  else if(is.vector(seqs)){
    seqs<-sapply(seqs,toupper)

    seqs_Lab<-alphabetCheck(seqs,alphabet = "dna",label)

    seqs<-seqs_Lab[[1]]
    label<-seqs_Lab[[2]]

  }
  else {
    stop("ERROR: Input sequence is not in the correct format. It should be a FASTA file or a string vector.")
  }


  lenSeqs<-sapply(seqs, nchar)

  numSeqs<-length(seqs)

  aaIdx<-c("A"=0.1260,"C"=0.1340,"G"=0.0806,"T"=0.1335)


featureMatrix<-kNUComposition_DNA(seqs=seqs,rng=3,reverse = FALSE,upto = FALSE,ORF = FALSE,normalized = FALSE)
  nam<-nameKmer(k=3,type = "dna")

  EIIP3mer<-vector(mode = "numeric",length = (4^3))

  for(j in 1:length(nam)){
    chars<-unlist(strsplit(nam[j],""))
    EIIP3mer[j]<-sum(aaIdx[chars])
  }

  for(i in 1:numSeqs){
    featureMatrix[i,]=featureMatrix[i,]*EIIP3mer
  }

  #featureMatrix<-t(t(featureMatrix)*EIIP3mer)

  if(length(label)==numSeqs){
    featureMatrix<-as.data.frame(featureMatrix)
    featureMatrix<-cbind(featureMatrix,label)
  }
  row.names(featureMatrix)<-names(seqs)



  return(featureMatrix)


}
