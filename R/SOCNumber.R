#' Sequence Order Coupling Number
#'
#' This function uses dissimilarity matrices Grantham and Schneider to compute the dissimilarity between amino acid pairs.
#' The distance between amino acid pairs is determined by d which varies between 1 to nlag.
#' For each d, it computes the sum of the dissimilarities of all amino acid pairs. The sum shows the value of tau for a value d.
#' The feature vector contains the values of taus for both matrices. Thus, the length of the feature vector is equal to nlag*2.
#'
#' @note When d=1, the pairs of amino acids have no gap and when d=2, there is one gap between the amino acid pairs in the sequence.
#' It will repeat likewise for other values of d.
#'
#' @param seqs is a FASTA file with amino acid sequences. Each sequence starts
#' with a '>' character. Also, seqs could be a string vector. Each element of the vector is a peptide/protein sequence.
#'
#'
#' @param nlag is a numeric value which shows the maximum distance between two amino acids.
#' Distances can be 1, 2, ..., or nlag. Defult is 30.
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @return It returns a feature matrix. The number of rows is equal to the number of sequences and
#' the number of columns is (nlag*2). For each distance d, there are two values. One value for Granthman and another one for Schneider distance.
#'
#'
#' @export
#'
#' @examples
#'
#' filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
#'
#' mat<-SOCNumber(seqs=filePrs,nlag=25)



SOCNumber<-function(seqs,nlag=30,label=c()){


  path.pack=system.file("extdata",package="ftrCOOL")
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
  numSeqs<-length(seqs)


  lenSeqs<-sapply(seqs, nchar)



  if (!all(lenSeqs>nlag)){
    stop("ERORR: 'nlag' should be less than the length of each sequence")
  }

  dismatrixSCHadr<-paste0(path.pack,"/Schneider-Wrede.txt")
  dismatrixGRAadr<-paste0(path.pack,"/Grantham.txt")

  dismatrixSCH<-read.delim(dismatrixSCHadr)
  dismatrixGRA<-read.delim(dismatrixGRAadr)

  row.names(dismatrixGRA)<-dismatrixGRA[,1]
  dismatrixGRA<-dismatrixGRA[,-1]
  row.names(dismatrixSCH)<-dismatrixSCH[,1]
  dismatrixSCH<-dismatrixSCH[,-1]
  dismatrixSCH<-as.matrix(dismatrixSCH)
  dismatrixGRA<-as.matrix(dismatrixGRA)
  featureMatrix<-matrix(0,nrow = numSeqs,ncol = (2*nlag))
  namesSCH<-paste("d=",1:nlag,"SCH")
  namesGRA<-paste("d=",1:nlag,"GRA")
  nameVect<-c(namesSCH,namesGRA)
  colnames(featureMatrix)<-nameVect
  for(n in 1:numSeqs){
    seq<-seqs[n]
    chars<-unlist(strsplit(seq,split = ""))
    N=nchar(seq)
    featureVector<-vector(mode="numeric",length = (2*nlag))
    vecti<-1:(N-1)
    vectj<-vecti+1
    TAUdSCHMAT<-dismatrixSCH[chars[vecti],chars[vectj]]
    TAUdGRAMAT<-dismatrixGRA[chars[vecti],chars[vectj]]
    for(d in nlag:1){
      TAUSCHd<-diag(TAUdSCHMAT)
      TAUSCHd<-TAUSCHd^2
      sumTAUSCHd<-sum(TAUSCHd)

      TAUdSCHMAT=TAUdSCHMAT[,-1]

      TAUGRAd<-diag(TAUdGRAMAT)
      TAUGRAd<-TAUGRAd^2
      sumTAUGRAd<-sum(TAUGRAd)

      TAUdGRAMAT=TAUdGRAMAT[,-1]

      featureVector[(nlag-d+1)]<-sumTAUSCHd
      featureVector[(2*nlag-d+1)]<-sumTAUGRAd

    }
    featureMatrix[n,]=featureVector
  }
  if(length(label)==numSeqs){
    featureMatrix<-as.data.frame(featureMatrix)
    featureMatrix<-cbind(featureMatrix,label)
  }
  row.names(featureMatrix)<-names(seqs)
  return(featureMatrix)

}
