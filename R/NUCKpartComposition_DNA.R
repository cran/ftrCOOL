#' Nucleotide to K Part Composition (NUCKpartComposition_DNA)
#'
#' In this function, each sequence is divided into k equal partitions.
#' The length of each part is equal to ceiling(l(lenght of the sequence)/k).
#' The last part can have a different length containing the residual nucleotides.
#' The nucleotide composition is calculated for each part.
#'
#' @note Warning: The length of all sequences should be greater than k.
#'
#' @param seqs is a FASTA file containing nucleotide sequences. The sequences start
#' with '>'. Also, seqs could be a string vector. Each element of the vector is a nucleotide sequence.
#'
#' @param k is an integer value. Each sequence should be divided to k partition(s).
#'
#' @param ORF (Open Reading Frame) is a logical parameter. If it is set to true, ORF region of each sequence is considered instead of the original sequence (i.e., 3-frame).
#'
#' @param reverseORF is a logical parameter. It is enabled only if ORF is true.
#' If reverseORF is true, ORF region will be searched in the sequence and also in the reverse complement of the sequence (i.e., 6-frame).
#'
#' @param normalized is a logical parameter. When it is FALSE, the return value of the function does not change. Otherwise, the return value is normalized using the length of the sequence.
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @return a feature matrix with k*4 number of columns. The number of rows is equal to the number of
#' sequences.
#'
#'
#' @export
#'
#' @examples
#'
#' fileLNC<-system.file("extdata/Athaliana_LNCRNA.fa",package="ftrCOOL")
#' mat<-NUCKpartComposition_DNA(seqs=fileLNC,k=5,ORF=TRUE,reverseORF=FALSE,normalized=FALSE)



NUCKpartComposition_DNA<- function(seqs,k=5,ORF=FALSE,reverseORF=TRUE,normalized=TRUE,label=c()) {


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

  if(ORF==TRUE){
    seqs=maxORF(seqs,reverse=reverseORF)
  }

  lenSeqs<-sapply(seqs,nchar)


  if(k<=0)
  {
    stop("k should be greater than zero")
  }

  if(!all(lenSeqs>k))
  {
    stop("ERROR: k should be smaller than length of sequences")
  }




  dict<-list("A"=1,"C"=2,"G"=3,"T"=4)


  numSeqs<-length(seqs)
  featureMatrix<-matrix(0,ncol = (4*k),nrow = numSeqs)

  tname<-nameKmer(k=1,type = "dna")
  tname<-rep(tname,k)
  wname<-rep(1:k,each =4)
  coln<-paste(tname,"p",wname,sep = "")

  colnames(featureMatrix)<-coln
  winSize<-ceiling(lenSeqs/k)

  for(n in 1:numSeqs){
    seq<-seqs[n]
    charSeq<-unlist(strsplit(seq,split = ""))

    for(i in 1:k)
    {

      winSeq<-charSeq[(((i-1)*winSize[n])+1):(i*winSize[n])]
      twinSeq<-table(winSeq)
      ntwinseq<-names(twinSeq)

      for(j in 1:length(twinSeq))
      {
        colindex<-(i-1)*4+as.numeric(dict[ntwinseq[j]])
        featureMatrix[n,colindex]<-twinSeq[j]

      }


    }
  }

  if(normalized==TRUE){
    featureMatrix[,]<-apply(featureMatrix, 2, function(i) i/winSize)

  }


  if(length(label)==numSeqs){
    featureMatrix<-as.data.frame(featureMatrix)
    featureMatrix<-cbind(featureMatrix,label)
  }
  row.names(featureMatrix)<-names(seqs)
  return(featureMatrix)

}
