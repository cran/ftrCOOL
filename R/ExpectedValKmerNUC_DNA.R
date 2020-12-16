#' Expected Value for K-mer Nucleotide (ExpectedValKmerNUC_DNA)
#'
#' This function is introduced by this package for the first time.
#' It computes the expected value for each k-mer in a sequence.
#' ExpectedValue(k-mer) = freq(k-mer) / ( freq(nucleotide1) * freq(nucleotide2)  * ... * freq(nucleotidek) )
#'
#' @param seqs is a FASTA file containing nucleotide sequences. The sequences start
#' with '>'. Also, seqs could be a string vector. Each element of the vector is a nucleotide sequence.
#'
#' @param k is an integer value. The default is four.
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
#' @return The function returns a feature matrix. The number of rows is equal to the number of sequences and
#' the number of columns is (4^k).
#'
#'
#' @export
#'
#' @examples
#'
#' fileLNC<-system.file("extdata/Athaliana_LNCRNA.fa",package="ftrCOOL")
#' mat<-ExpectedValKmerNUC_DNA(seqs=fileLNC,k=4,ORF=TRUE,reverseORF=FALSE)



ExpectedValKmerNUC_DNA<-function(seqs,k=4,ORF=FALSE,reverseORF=TRUE,normalized=TRUE,label=c())
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
  flag=0
  if(ORF==TRUE){
    if(length(label)==length(seqs)){
      names(label)=names(seqs)
      flag=1
    }
    seqs=maxORF(seqs,reverse=reverseORF)
    if(flag==1)
      label=label[names(seqs)]
  }
  numSeqs<-length(seqs)
  featureMatrix<-kNUComposition_DNA(seqs,reverse=FALSE,rng=k,normalized=FALSE)
  NUCompos<-kNUComposition_DNA(seqs,rng=1,reverse=FALSE,normalized=FALSE)

  for(i in 1:ncol(featureMatrix)){
    chars<-unlist(strsplit(colnames(featureMatrix)[i],split = ""))
    composMat<-NUCompos[,chars]
    if(is.matrix(composMat)){
      mult<-apply(composMat, 1, prod)
    } else {
      mult<-prod(composMat)
    }
    featureMatrix[,i]<-featureMatrix[,i]/mult
  }
  featureMatrix[is.na(featureMatrix)]=0
  featureMatrix[is.infinite(featureMatrix)]=0

  if(normalized==TRUE){
    seqLen<-sapply(seqs, nchar)
    featureMatrix<-featureMatrix/seqLen
  }

  if(length(label)==numSeqs){
    featureMatrix<-as.data.frame(featureMatrix)
    featureMatrix<-cbind(featureMatrix,label)
  }
  row.names(featureMatrix)<-names(seqs)
  return(featureMatrix)
}
