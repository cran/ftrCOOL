#' Expected Value for each Amino Acid (ExpectedValueAA)
#'
#' This function is introduced by this package for the first time.
#' It computes the expected value for each k-mer in a sequence.
#' ExpectedValue(k-mer) = freq(k-mer) / (c_1 * c_2  * ... * c_k),
#' where c_i is the number of codons that encrypt the i'th amino acid in the k-mer.
#'
#' @param seqs is a FASTA file with amino acid sequences. Each sequence starts
#' with a '>' character. Also, seqs could be a string vector. Each element of the vector is a peptide/protein sequence.
#'
#' @param k is an integer value. The default is two.
#'
#' @param normalized is a logical parameter. When it is FALSE, the return value of the function does not change. Otherwise, the return value is normalized using the length of the sequence.
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @return This function returns a feature matrix. The number of rows is equal to the number of sequences and
#' the number of columns is 20^k.
#'
#'
#'
#' @export
#'
#' @examples
#'
#' filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
#' mat<-ExpectedValueAA(seqs=filePrs,k=2,normalized=FALSE)





ExpectedValueAA<-function(seqs,k=2,normalized=TRUE,label=c())
{

  if(length(seqs)==1&&file.exists(seqs)){
    seqs<-fa.read(seqs,alphabet="aa")
    seqs_Lab<-alphabetCheck(seqs,alphabet = "aa",label)

    seqs<-seqs_Lab[[1]]
    label<-seqs_Lab[[2]]
  }
  else if(is.vector(seqs)){
    seqs<-sapply(seqs,toupper)

    seqs_Lab<-alphabetCheck(seqs,alphabet = "aa",label)

    seqs<-seqs_Lab[[1]]
    label<-seqs_Lab[[2]]

  }else {
    stop("ERROR: Input sequence is not in the correct format. It should be a FASTA file or a string vector.")
  }
  numSeqs<-length(seqs)
  featureMatrix<-kAAComposition(seqs,rng=k,normalized = FALSE)

  numCodAA <- list("A"=4,"C"=2,"D"=2,"E"=2,"F"=2,"G"=4,"H"=2,"I"=3,"K"=2,"L"=6,"M"=1,"N"=2,"P"=4,"Q"=2,"R"=6,"S"=6,"T"=4,"V"=4,"W"=1,"Y"=2)




  for(i in 1:ncol(featureMatrix)){
    chars<-unlist(strsplit(colnames(featureMatrix)[i],split = ""))
    numAA<-unlist(numCodAA[chars])
    mult<-prod(numAA)
    featureMatrix[,i]<-featureMatrix[,i]/mult
  }

  featureMatrix[is.na(featureMatrix)]=0

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
