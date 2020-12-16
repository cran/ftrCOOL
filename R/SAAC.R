#' Splitted Amino Acid Composition (SAAC)
#'
#' This function splits the input sequence into three parts. The first part is N-terminal
#' and the third part is C-terminal and middle part contains all amino acids between these two part.
#' N-terminal will be determined by the first numNterm amino acid in the sequences and C-terminal is
#' determined by numCterm of the last amino acids in the sequence. Users should enter numNterm and numCterm
#' parameters. Their default value is 25.
#' The function calculates \link{kAAComposition} for each of the three parts.
#'
#'
#' @param seqs is a FASTA file with amino acid sequences. Each sequence starts
#' with a '>' character. Also, seqs could be a string vector. Each element of the vector is a peptide/protein sequence.
#'
#'
#' @param k shows which type of amino acid composition applies to the parts.
#' For example, the amino acid composition is applied when k=1 and when k=2, the dipeptide Composition is applied.
#'
#' @param numNterm shows how many amino acids should be considered for N-terminal.
#'
#' @param numCterm shows how many amino acids should be considered for C-terminal.
#'
#' @param normalized is a logical parameter. When it is FALSE, the return value of the function does not change. Otherwise, the return value is normalized using the length of the sequence.
#'
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @return It returns a feature matrix. The number of rows is equal to the number of sequences.
#' The number of columns is (3*(20^k)).
#'
#'
#'
#' @export
#'
#' @examples
#'
#' filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
#' mat<-SAAC(seqs=filePrs,k=1,numNterm=15,numCterm=15)


SAAC <- function(seqs,k=1,numNterm=5,numCterm=5,normalized=TRUE,label=c()) {

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




  lenSeqs<-sapply(seqs, nchar)
  if(!all(lenSeqs>(numNterm+numCterm))){
    deletInd<-which(lenSeqs<(numNterm+numCterm))
    deletedNames<-names(deletInd)
    strNames<-toString(deletedNames)
    lens<-lenSeqs[deletInd]
    strlens<-toString(lens)

    if(length(deletInd)==length(seqs)){
      stop("ERROR: All sequences were deleted!")
    }

    warning(paste("Sequences",strNames,"with lengths",strlens,"were deleted. The aggregate lenght of Nterminal and Cterminal parts in a sequence should be less than the length of the sequence."))

    if(length(label)==length(lenSeqs)){
      label<-label[-deletInd]
    }
    lenSeqs<-lenSeqs[-deletInd]
    seqs<-seqs[-deletInd]

  }
  numSeqs<-length(seqs)

  ntermSeq<-substr(seqs,1,numNterm)
  ctermSeq<-substr(seqs,(lenSeqs-numCterm+1),lenSeqs)
  midtermSeq<-substr(seqs,(numNterm+1),(lenSeqs-numCterm))
  nTermMat<-kAAComposition(ntermSeq,rng=k,normalized = FALSE)
  cTermMat<-kAAComposition(ctermSeq,rng=k,normalized = FALSE)
  midTermMat<-kAAComposition(midtermSeq,rng=k,normalized = FALSE)
  if(normalized==TRUE){
    nTermMat<-nTermMat/numNterm
    cTermMat<-cTermMat/numCterm
    lenMidle<-(lenSeqs-(numNterm+numCterm))
    midTermMat<-midTermMat/lenMidle
  }
  colnames(nTermMat)<-paste("N_",colnames(nTermMat),sep = "")
  colnames(cTermMat)<-paste("C_",colnames(cTermMat),sep = "")
  colnames(midTermMat)<-paste("M_",colnames(midTermMat),sep = "")


  featureMatrix<-cbind(nTermMat,midTermMat,cTermMat)
  if(length(label)==numSeqs){
    featureMatrix<-as.data.frame(featureMatrix)
    featureMatrix<-cbind(featureMatrix,label)
  }
  row.names(featureMatrix)<-names(seqs)
  return(featureMatrix)


}
