#' Amino Acid to K Part Composition (AAKpartComposition)
#'
#' In this function, each sequence is divided into k equal partitions.
#' The length of each part is equal to ceiling(l(lenght of the sequence)/k).
#' The last part can have a different length containing the residual amino acids.
#' The amino acid composition is calculated for each part.
#'
#' @note Warning: The length of all sequences should be greater than k.
#'
#'
#' @param seqs is a FASTA file with amino acid sequences. Each sequence starts
#' with a '>' character. Also, seqs could be a string vector. Each element of the vector is a peptide/protein sequence.
#'
#' @param k is an integer value. Each sequence should be divided to k partition(s).
#'
#' @param normalized is a logical parameter. When it is FALSE, the return value of the function does not change. Otherwise, the return value is normalized using the length of the sequence.
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @return a feature matrix with k*20 number of columns. The number of rows is equal to the number of
#' sequences.
#'
#' @export
#'
#' @examples
#'
#' filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
#' mat<-AAKpartComposition(seqs=filePrs,k=5,normalized=FALSE)
#'

AAKpartComposition<- function(seqs,k=3,normalized=TRUE,label=c()) {


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

  }
  else {
    stop("ERROR: Input sequence is not in the correct format. It should be a FASTA file or a string vector.")
  }

  lenSeq<-sapply(seqs,nchar)

  if(k<=0)
  {
    stop("ERROR: k should be greater than zero")
  }



  if(!all(lenSeq>=k)){
    deletInd<-which(lenSeq<k)
    deletedNames<-names(deletInd)
    strNames<-toString(deletedNames)
    lens<-lenSeq[deletInd]
    strlens<-toString(lens)

    warning(paste("Sequences",strNames,"with lengths",strlens,"were deleted. Their lenghts were less than k"))

    if(length(label)==length(lenSeq)){
      label<-label[-deletInd]
    }
    lenSeq<-lenSeq[-deletInd]
    seqs<-seqs[-deletInd]
  }



  dict<-list("A"=1,"C"=2,"D"=3,"E"=4,"F"=5,"G"=6,"H"=7,"I"=8,"K"=9,"L"=10,"M"=11,"N"=12,"P"=13,"Q"=14,"R"=15,"S"=16,"T"=17,"V"=18,"W"=19,"Y"=20)


  numSeqs<-length(seqs)
  featureMatrix<-matrix(0,ncol = (20*k),nrow = numSeqs)

  tname<-nameKmer(k=1,type = "aa")
  tname<-rep(tname,k)
  wname<-rep(1:k,each =20)
  coln<-paste(tname,"p",wname,sep = "")

  colnames(featureMatrix)<-coln
  winSize<-ceiling(lenSeq/k)

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
        colindex<-(i-1)*20+as.numeric(dict[ntwinseq[j]])
        featureMatrix[n,colindex]<-twinSeq[j]
      }


    }
  }

  if(normalized==TRUE){
    featureMatrix[,]<-apply(featureMatrix, 2, function(i) i/lenSeq)
  }
  row.names(featureMatrix)<-names(seqs)

  if(length(label)==numSeqs){
    featureMatrix<-as.data.frame(featureMatrix)
    featureMatrix<-cbind(featureMatrix,label)
  }

  return(featureMatrix)

}



