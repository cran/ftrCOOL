#' Dipeptide Deviation from Expected Mean value
#'
#' This function computes the dipeptide deviation from the expected mean value.
#'
#'
#'
#' @param seqs is a FASTA file with amino acid sequences. Each sequence starts
#' with a '>' character. Also, seqs could be a string vector. Each element of the vector is a peptide/protein sequence.
#'
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @return A feature matrix with 20^2=400 number of columns. The number of rows is equal to the number of sequences.
#'
#'
#' @export
#'
#' @examples
#'
#' filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
#' mat<-DDE(seqs=filePrs)
#'



DDE <-function(seqs,label=c())
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

  #is a matrix with ncol=400
  Dc <-Dc(seqs)

  # is a vector with length 400
  Tm <-Tm()
  names(Tm)<-colnames(Dc)

  charSeqs<-sapply(seqs,nchar)

  #is matrix with ncol=400
  Tv <-Tv(Tm,charSeqs)

  dde<-t(apply(Dc, 1, function(i) i-Tm))

  #is a matrix with ncol=(lenGrp)^2
  dde<-dde/sqrt(Tv)


  colnames(dde)<-nameKmer(k=2,type = "aa")

  if(length(label)==numSeqs){
    dde<-as.data.frame(dde)
    dde<-cbind(dde,label)
  }
  row.names(dde)<-names(seqs)

  return(dde)

}

Dc <- function(seqs)
{
  len<-sapply(seqs, nchar)
  len<-len-1
  dipepCompos <- kAAComposition(seqs,rng=2,upto=FALSE,normalized = FALSE)
  dipepCompos<-dipepCompos/len
  return(dipepCompos)
}

Tm <- function(){

  numCodAA <- list("A"=4,"C"=2,"D"=2,"E"=2,"F"=2,"G"=4,"H"=2,"I"=3,"K"=2,"L"=6,"M"=1,"N"=2,"P"=4,"Q"=2,"R"=6,"S"=6,"T"=4,"V"=4,"W"=1,"Y"=2)


  Cn<- sum(unlist(numCodAA))
  tm <- matrix(nrow =20 ,ncol = 20)

  for (i in 1:20){
    for(j in 1:20){
      tm[j,i]<-as.numeric(numCodAA[i])*as.numeric(numCodAA[j])/(Cn*Cn)
    }

  }
  tm<-as.vector(tm)

  return(tm)
}

Tv <- function(tm,lens)
{
  x<-(tm*(1-tm))
  TvMatrix<-matrix(0,ncol = length(x),nrow = length(lens))
  colnames(TvMatrix)<-colnames(x)
  for(i in 1:length(lens)){
    TvMatrix[i,]<-x/(lens[i]-1)
  }
  return(TvMatrix)

}

