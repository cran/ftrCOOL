#' Expected Value for K-mer Amino Acid (ExpectedValueKmerAA)
#'
#' This function computes the expected value of each k-mer by dividing the frequency of the
#' kmer to multiplying frequency of each amino acid of the k-mer in the sequence.
#'
#' ExpectedValue(k-mer) = freq(k-mer) / ( freq(aminoacid1) * freq(aminoacid2)  * ... * freq(aminoacidk) )
#'
#' @param seqs is a FASTA file with amino acid sequences. Each sequence starts
#' with a '>' character. Also, seqs could be a string vector. Each element of the vector is a peptide/protein sequence.
#'
#' @param k is an integer value and it shows the size of kmer in the kmer composition. The default value is 2.
#'
#' @param normalized is a logical parameter. When it is FALSE, the return value of the function does not change. Otherwise, the return value is normalized using the length of the sequence.
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @return This function returns a feature matrix. The number of rows equals the number of sequences and
#' the number of columns if upto set false, is 20^k.
#'
#'
#' @export
#'
#' @examples
#'
#' filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
#' mat<-ExpectedValueKmerAA(filePrs,k=2,normalized=FALSE)


ExpectedValueKmerAA<-function(seqs,k=2,normalized=TRUE,label=c())
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
  featureMatrix<-kAAComposition(seqs,rng=k,normalized =FALSE)
  AAcompos<-kAAComposition(seqs,rng=1,normalized =FALSE)


  for(i in 1:ncol(featureMatrix)){

    chars<-unlist(strsplit(colnames(featureMatrix)[i],split = ""))
    composMat<-AAcompos[,chars]
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
