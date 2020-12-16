#' Composition of k-Spaced Amino Acids pairs (CkSAApair)
#'
#'
#' This function calculates the composition of k-spaced amino acid pairs. In other words, it
#' computes the frequency of all amino acid pairs with k spaces.
#'
#' @note 'upto' is enabled only when rng is a number and not a vector.
#'
#' @param seqs is a FASTA file with amino acid sequences. Each sequence starts
#' with a '>' character. Also, seqs could be a string vector. Each element of the vector is a peptide/protein sequence.
#'
#' @param rng This parameter can be a number or a vector. Each element of the vector shows the number of spaces between amino acid pairs.
#' For each k in the rng vector, a new vector (whose size is 400) is created which contains the frequency of pairs with k gaps.
#'
#' @param upto It is a logical parameter. The default value is FALSE. If rng is a number and upto is set to TRUE, rng is converted
#' to a vector with values from [0 to rng].
#'
#' @param normalized is a logical parameter. When it is FALSE, the return value of the function does not change. Otherwise, the return value is normalized using the length of the sequence.
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @return The function returns a feature matrix. The number of rows is equal to the number of sequences and
#' the number of columns is 400*(length of rng vector).
#'
#'
#' @export
#'
#' @examples
#'
#' filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
#' mat1<-CkSAApair(seqs=filePrs,rng=2,upto=TRUE,normalized=TRUE)
#'
#' mat2<-CkSAApair(seqs=filePrs,rng=c(1,3,5))
#'

CkSAApair <- function(seqs,rng=3,upto=FALSE,normalized=TRUE,label=c()){


  dict<-list("A"=1,"C"=2,"D"=3,"E"=4,"F"=5,"G"=6,"H"=7,"I"=8,"K"=9,"L"=10,"M"=11,"N"=12,"P"=13,"Q"=14,"R"=15,"S"=16,"T"=17,"V"=18,"W"=19,"Y"=20)


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

  numSeqs=length(seqs)

  if(upto==TRUE && length(rng)==1){
    l<-length(rng)
    l<-rng[l]
    rng<-0:l
  }


  rng <- sort(rng)
  rng <- unique(rng)
  len<-length(rng)

  featureMatrix <- matrix(0 , ncol = len*400,nrow = numSeqs)
  dipep<-nameKmer(k=2,type = "aa")
  for(i in 1:length(dipep)){
    ditemp<-unlist(strsplit(dipep[i],split = ""))
    dipep[i]<-paste(ditemp[1],ditemp[2])
  }


  featName<-vector()
  for(i in 1:len){
    featName<-c(featName,gsub(" ",strrep("s",rng[i]),dipep))
  }

  colnames(featureMatrix)<-featName



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
      for(j in 1:length(tbkmers)){
        tmp<-unlist(strsplit(nmtbkmers[j],split = ""))
        index<-(as.numeric(dict[tmp[1]])-1)*20+as.numeric(dict[tmp[2]])
        index<-index+((i-1)*400)
        featureMatrix[n,index]<-tbkmers[j]
      }
    }

  }

  if(normalized==TRUE){
    seqLen<-sapply(seqs, nchar)
    featureMatrix<-featureMatrix/seqLen
  }
  row.names(featureMatrix)<-names(seqs)
  if(length(label)==numSeqs){
    featureMatrix<-as.data.frame(featureMatrix)
    featureMatrix<-cbind(featureMatrix,label)
  }

  return(featureMatrix)
}



