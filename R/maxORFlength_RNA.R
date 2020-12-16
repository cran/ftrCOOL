#' Maximum Open Reading Frame length in RNA (maxORFlength_RNA)
#'
#' This function returns the length of the maximum Open Reading Frame for each sequence.
#' If reverse is FALSE, ORF region will be searched in a sequence.
#' Otherwise, it will be searched both in the sequence and its reverse complement.
#'
#'
#' @param seqs is a FASTA file containing ribonucleotide sequences. The sequences start
#' with '>'. Also, seqs could be a string vector. Each element of the vector is a ribonucleotide sequence.
#'
#' @param reverse It is a logical parameter which assumes the reverse complement of the sequence.
#'
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @param normalized is a logical parameter. When it is FALSE, the return value of the function does not change. Otherwise, the return value is normalized using the length of the sequence.
#'
#'
#' @return A vector containing the lengths of maximum ORFs for each sequence.
#'
#'
#'
#' @export
#'
#' @examples
#'
#' fileLNC<-system.file("extdata/Carica_papaya101RNA.txt",package="ftrCOOL")
#' vect<-maxORFlength_RNA(seqs=fileLNC,reverse=TRUE,normalized=TRUE)
#'

maxORFlength_RNA<-function(seqs,reverse=TRUE,normalized=FALSE,label=c()){

  if(length(seqs)==1&&file.exists(seqs)){
    seqs<-fa.read(seqs,alphabet="rna")
    seqs_Lab<-alphabetCheck(seqs,alphabet = "rna",label)

    seqs<-seqs_Lab[[1]]
    label<-seqs_Lab[[2]]
  }
  else if(is.vector(seqs)){
    seqs<-sapply(seqs,toupper)

    seqs_Lab<-alphabetCheck(seqs,alphabet = "rna",label)

    seqs<-seqs_Lab[[1]]
    label<-seqs_Lab[[2]]

  }
  else {
    stop("ERROR: Input sequence is not in the correct format. It should be a FASTA file or a string vector.")
  }
  numSeqs<-length(seqs)
  listORF<-list()
  listSeqs<-list()
  orfLen<-vector(mode = "numeric",length = numSeqs)
  for(n in 1:numSeqs){
    seq<-seqs[n]
    len<-nchar(seq)
    firstSeq<-substring(seq,seq(1,(len-2),3),seq(3,len,3))
    secSeq<-substring(seq,seq(2,(len-3),3),seq(4,len,3))
    thSeq<-substring(seq,seq(3,(len-4),3),seq(5,len,3))
    listSeqs[[1]]<-firstSeq
    listSeqs[[2]]<-secSeq
    listSeqs[[3]]<-thSeq
    listORF[[1]]<-findORF_RNA(firstSeq)
    listORF[[2]]<-findORF_RNA(secSeq)
    listORF[[3]]<-findORF_RNA(thSeq)

    if(reverse==TRUE){
      revSeq<-revComp(seq)
      revSeq<-paste(revSeq,collapse = "")
      firstSeq<-substring(revSeq,seq(1,(len-2),3),seq(3,len,3))
      secSeq<-substring(revSeq,seq(2,(len-3),3),seq(4,len,3))
      thSeq<-substring(revSeq,seq(3,(len-4),3),seq(5,len,3))
      listORF[[4]]<-findORF_RNA(firstSeq)
      listORF[[5]]<-findORF_RNA(secSeq)
      listORF[[6]]<-findORF_RNA(thSeq)

    }
    orfVect=c()
    for(i in 1:length(listORF)){
      orfVect<-c(orfVect,listORF[[i]][1])
    }

    maxORF<-max(orfVect)

    orfLen[n]<-(maxORF*3)


  }
  if(normalized==TRUE){
    lenSeqs<-sapply(seqs, nchar)
    orfLen<-orfLen/lenSeqs
  }
  names(orfLen)<-names(seqs)
  if(length(label)==numSeqs){

    orfLen<-cbind(orfLen,label)
  }
  return(orfLen)
}

# findORF_RNA<-function(Triplet){
#   strCds=match("AUG",Triplet)
#   stopCds=0
#   lenORF=0
#   orf<-c(0,0,0)
#   names(orf)<-c("lenORF","indStr","indStp")
#   if(!is.na(strCds)){
#     mtc=match(c("UAG","UAA","UGA"),Triplet)
#     mtc[is.na(mtc)]<-0
#     inds<-which(mtc>strCds)
#     if(length(inds)>0){
#       stopCds<-min(mtc[inds])
#       lenORF=stopCds-strCds+1
#       orf<-c(lenORF,strCds,stopCds)
#     }
#     else{
#       lenORF<-0
#     }
#
#   }
#   return(orf)
#
# }
