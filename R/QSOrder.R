#' Quasi Sequence Order
#'
#' This function computes the quasi-sequence-order for sequences.
#' It is for amino acid pairs with d distances (d can be any number between 1 and 20).
#' First, it calculates the frequencies of each amino acid ("A", "C",..., "Y").
#' Then, it normalizes the frequencies by dividing the frequency of an amino acid to the frequency of all amino acids
#' plus the sum of tau values which is multiplied by W. tau values are given by function \link{SOCNumber}.
#' For d bigger than 20, it computes tau for d in the range "1 to (nlag-20) * W" and normalizes them like before.
#'
#' @details Please find details about tau in function \link{SOCNumber}.
#'
#' @note For d between 21 to nlag, the function calculates tau values for (d-20) to (nlag-20).
#'
#' @param seqs is a FASTA file with amino acid sequences. Each sequence starts
#' with a '>' character. Also, seqs could be a string vector. Each element of the vector is a peptide/protein sequence.
#'
#' @param nlag is a numeric value which shows the maximum distance between two amino acids.
#' Distances can be 1, 2, ..., or nlag.
#'
#' @param W (weight) is a tuning parameter.
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @return It returns a feature matrix which the number of rows equals to the number of sequences and the number of columns is (nlag*2).
#' For each distance d, there are two values. One value for Granthman and another one for Schneider distance.
#'
#'
#' @export
#'
#' @examples
#'
#' filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
#'
#' mat<-QSOrder(seqs=filePrs,nlag=25)
#'
QSOrder<-function(seqs,nlag=25,W=0.1,label=c()){


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

  seqs<-sapply(seqs, toupper)
  lenSeqs<-sapply(seqs, nchar)

  if (!all(lenSeqs>nlag)){
    stop("ERORR: 'nlag' should be less than the length of each sequence")
  }



  TAU<-SOCNumber(seqs,nlag)
  nTAU<-nrow(TAU)
  TauSCH<-matrix(ncol = nlag,nrow = nTAU)
  TauGRA<-matrix(ncol = nlag,nrow = nTAU)
  TauSCH[,1:nlag]<-TAU[,1:nlag]
  TauGRA[,1:nlag]<-TAU[,(nlag+1):(2*nlag)]
  f<-kAAComposition(seqs = seqs,rng=1,normalized = FALSE)


    sigmaTauSCH<-apply(TauSCH, 1, sum)
    sigmaTauGRA<-apply(TauGRA, 1, sum)
    sigmafr<-apply(f, 1, sum)




  xrSCH<-matrix(nrow = nTAU,ncol = 20)
  xrGRA<-matrix(nrow = nTAU,ncol = 20)

  xrSCH<-f/(sigmafr+W*sigmaTauSCH)
  xrGRA<-f/(sigmafr+W*sigmaTauGRA)


  colnames(xrSCH)<-paste("SCH",paste("f",colnames(f),sep = ""))
  colnames(xrGRA)<-paste("GRA",paste("f",colnames(f),sep = ""))
  xr2SCH=c()
  xr2GRA=c()
  if(nlag>20){
    xr2SCH=matrix(nrow = nTAU,ncol = (nlag-20))
    xr2GRA=matrix(nrow = nTAU,ncol = (nlag-20))
    xr2SCH[,]=((W*TauSCH[,21:nlag])-20)/(sigmafr+W*sigmaTauSCH)
    xr2GRA[,]=((W*TauGRA[,21:nlag])-20)/(sigmafr+W*sigmaTauGRA)
    colnames(xr2SCH)=paste("SCH d=",21:nlag)
    colnames(xr2GRA)=paste("GRA d=",21:nlag)
  }
  featureMatrix<-cbind(xrSCH,xr2SCH,xrGRA,xr2GRA)
  if(length(label)==numSeqs){
    featureMatrix<-as.data.frame(featureMatrix)
    featureMatrix<-cbind(featureMatrix,label)
  }
  row.names(featureMatrix)<-names(seqs)
  return(featureMatrix)

}
