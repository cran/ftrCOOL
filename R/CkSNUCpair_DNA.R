#' Composition of k-Spaced Nucleotides Pairs (CkSNUCpair_DNA)
#'
#' This function calculates the composition of k-spaced nucleotide pairs. In other words, it
#' computes the frequency of all nucleotide pairs with k spaces.
#'
#' @note 'upto' is enabled only when rng is a number and not a vector.
#'
#' @param seqs is a FASTA file containing nucleotide sequences. The sequences start
#' with '>'. Also, seqs could be a string vector. Each element of the vector is a nucleotide sequence.
#'
#' @param rng This parameter can be a number or a vector. Each element of the vector shows the number of spaces between nucleotide pairs.
#' For each k in the rng vector, a new vector (whose size is 16) is created which contains the frequency of pairs with k gaps.
#'
#' @param upto It is a logical parameter. The default value is FALSE. If rng is a number and upto is set to TRUE, rng is converted
#' to a vector with values from [0 to rng].
#'
#'
#' @param ORF (Open Reading Frame) is a logical parameter. If it is set to true, ORF region of each sequence is considered instead of the original sequence (i.e., 3-frame).
#'
#' @param reverseORF is a logical parameter. It is enabled only if ORF is true.
#' If reverseORF is true, ORF region will be searched in the sequence and also in the reverse complement of the sequence (i.e., 6-frame).
#'
#'
#' @param normalized is a logical parameter. When it is FALSE, the return value of the function does not change. Otherwise, the return value is normalized using the length of the sequence.
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @return The function returns a feature matrix. The number of rows is equal to the number of sequences and
#' the number of columns is 16*(length of rng vector).
#'
#'
#' @export
#'
#' @examples
#'
#'
#' fileLNC<-system.file("extdata/Athaliana_LNCRNA.fa",package="ftrCOOL")
#' mat1<-CkSNUCpair_DNA(seqs=fileLNC,rng=2,upto=TRUE,ORF=TRUE,reverseORF=FALSE)
#' mat2<-CkSNUCpair_DNA(seqs=fileLNC,rng=c(1,3,5))
#'



CkSNUCpair_DNA <- function(seqs,rng=3,upto=FALSE,ORF=FALSE,reverseORF=TRUE,normalized=TRUE,label=c()){



  dict<-list("A"=1,"C"=2,"G"=3,"T"=4)


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

  numSeqs=length(seqs)

  if(upto==TRUE && length(rng)==1){
    l<-length(rng)
    l<-rng[l]
    rng<-0:l
  }


  rng <- sort(rng)
  len<-length(rng)

  featureMatrix <- matrix(0 , ncol = len*16,nrow = numSeqs)
  dipep<-nameKmer(k=2,type = "dna")
  for(i in 1:length(dipep)){
    ditemp<-unlist(strsplit(dipep[i],split = ""))
    dipep[i]<-paste(ditemp[1],ditemp[2])
  }


  featName<-vector()
  for(i in 1:len){
    featName<-c(featName,gsub(" ",strrep("s",rng[i]),dipep))
  }

  colnames(featureMatrix)<-featName
  tempname<-nameKmer(k=2,type = "dna")


  for(n in 1:numSeqs){
    seq<-seqs[n]
    seqChars<-unlist(strsplit(seq,split = ""))
    lenSeq<-length(seqChars)

    for(i in 1:len){
      temp1<-seqChars[1:(lenSeq-rng[i]-1)]
      temp2<-seqChars[((rng[i]+1)+1):(lenSeq)]
      kmers<-paste(temp1,temp2,sep = "")
      tbkmers<-table(kmers)
      nmtbkmers<-names(tbkmers)


      ((len-1)*16)+1:len*16
      tempvect<-vector(mode = "numeric",length = 16)
      names(tempvect)<-tempname
      tempvect[nmtbkmers]<-tbkmers
      featureMatrix[n,(((len-1)*16)+1):(len*16)]<-tempvect
    }

  }

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



