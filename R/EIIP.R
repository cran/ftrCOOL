#' Electron-ion interaction pseudopotentials (EIIP)
#'
#' This function replaces each nucleotide in the input sequence with its electron-ion interaction value.
#' The resulting sequence is represented by a feature vector whose length is equal to the length of the sequence.
#' Please check the references for more information.
#'
#' @references Chen, Zhen, et al. "iLearn: an integrated platform and meta-learner for feature engineering, machine-learning analysis and modeling of DNA, RNA and protein sequence data." Briefings in bioinformatics 21.3 (2020): 1047-1057.
#'
#' @param seqs is a FASTA file containing nucleotide sequences. The sequences start
#' with '>'. Also, seqs could be a string vector. Each element of the vector is a nucleotide sequence.
#'
#'
#' @param outFormat (output format) can take two values: 'mat'(matrix) and 'txt'. The default value is 'mat'.
#'
#' @param outputFileDist shows the path and name of the 'txt' output file.
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @return The output depends on the outFormat parameter which can be either 'mat' or 'txt'. If outFormat is 'mat', the function returns a feature
#' matrix for sequences with the same length such that the number of columns is equal to the length of the sequences
#' and the number of rows is equal to the number of sequences. It is usable for machine learning purposes.
#' If the outFormat is 'txt', the output is written to a tab-delimited file.
#'
#'
#' @export
#'
#' @examples
#'
#' LNCSeqsADR<-system.file("extdata/",package="ftrCOOL")
#' LNC50Nuc<-as.vector(read.csv(paste0(LNCSeqsADR,"/LNC50Nuc.csv"))[,2])
#' mat<-EIIP(seqs = LNC50Nuc,outFormat="mat")


EIIP <- function(seqs,outFormat="mat",outputFileDist="",label=c())
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


  ###deleteing high correlate features

if(outFormat=="mat"){
  if(length(unique(lenSeqs))>1){
    stop("ERROR: All sequences should have the same length in 'mat' mode. For sequences with different lengths, please use 'txt' for outFormat parameter")
  }


  featureMatrix=matrix(0,nrow = numSeqs,ncol = lenSeqs[1])
  colna<-paste0("pos",1:lenSeqs[1])
  colnames(featureMatrix)<-colna

  for(n in 1:numSeqs){
    seq<-seqs[n]
    chars<-unlist(strsplit(seq,NULL))
    aaMatselected<-aaIdx[chars]
    featureMatrix[n,]<-aaMatselected
  }


  if(length(label)==numSeqs){
    featureMatrix<-as.data.frame(featureMatrix)
    featureMatrix<-cbind(featureMatrix,label)
  }
  row.names(featureMatrix)<-names(seqs)
  return(featureMatrix)
}
  else{

    nameSeq<-names(seqs)

    for(n in 1:numSeqs){
      seq<-seqs[n]
      chars<-unlist(strsplit(seq,NULL))
      aaMatselected<-aaIdx[chars]
      temp<-c(nameSeq[n],aaMatselected)
      temp2<-paste(temp,collapse = "\t")
      write(temp2,outputFileDist,append = TRUE)
    }



  }
}
