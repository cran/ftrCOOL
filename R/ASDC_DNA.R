#' Adaptive skip dinucleotide composition_DNA) (ASDC_DNA)
#'
#' This descriptor sufficiently considers the correlation information present not only between adjacent nucleotides but also between intervening nucleotides
#' This function calculates frequency of pair nucleotides omitting gaps between them. Then this function normalizes each value through dividing each frequency by summition(frequencies).
#'
#' @references Wei L, Zhou C, Chen H, Song J, Su R. ACPred-FL: a sequence-based predictor using effective feature representation to improve the prediction of anti-cancer peptides. Bioinformatics (2018).
#'
#' @param seqs is a FASTA file containing nucleotide sequences. The sequences start
#' with '>'. Also, seqs could be a string vector. Each element of the vector is a nucleotide sequence.
#'
#'
#'
#' @param ORF (Open Reading Frame) is a logical parameter. If it is set to true, ORF region of each sequence is considered instead of the original sequence (i.e., 3-frame).
#'
#' @param reverseORF is a logical parameter. It is enabled only if ORF is true.
#' If reverseORF is true, ORF region will be searched in the sequence and also in the reverse complement of the sequence (i.e., 6-frame).
#'
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @return The function returns a feature matrix. The number of rows is equal to the number of sequences and
#' the number of columns is 16 (All posible nucleotide pairs).
#'
#'
#' @export
#'
#' @examples
#'
#'
#' fileLNC<-system.file("extdata/Athaliana_LNCRNA.fa",package="ftrCOOL")
#' fileLNC<-fa.read(file=fileLNC,alphabet="dna")[1:5]
#' mat1<-ASDC_DNA(seqs=fileLNC,ORF=TRUE,reverseORF=FALSE)
#'



ASDC_DNA <- function(seqs,ORF=FALSE,reverseORF=TRUE,label=c()){


  upto=TRUE
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

  rng<-sapply(seqs, nchar)
  rng<-rng-1



  featureMatrix <- matrix(0 , ncol = 16,nrow = numSeqs)
  dipep<-nameKmer(k=2,type = "dna")
  # for(i in 1:length(dipep)){
  #   ditemp<-unlist(strsplit(dipep[i],split = ""))
  #   dipep[i]<-paste(ditemp[1],ditemp[2])
  # }


  featName<-vector()
  # for(i in 1:len){
  #   featName<-c(featName,gsub(" ",strrep("s",rng[i]),dipep))
  # }

  colnames(featureMatrix)<-dipep
  tempname<-dipep



  for(n in 1:numSeqs){
    seq<-seqs[n]
    seqChars<-unlist(strsplit(seq,split = ""))
    lenSeq<-length(seqChars)
    len<-rng[n]

    for(i in 1:len){
      nums<-i-1
      temp1<-seqChars[1:(lenSeq-nums-1)]
      temp2<-seqChars[((nums+1)+1):(lenSeq)]
      kmers<-paste(temp1,temp2,sep = "")
      tbkmers<-table(kmers)
      nmtbkmers<-names(tbkmers)


      tempvect<-vector(mode = "numeric",length = 16)
      names(tempvect)<-tempname
      tempvect[nmtbkmers]<-tbkmers
      featureMatrix[n,]<-featureMatrix[n,]+tempvect
    }
    featureMatrix[n,]<-featureMatrix[n,]/sum(featureMatrix[n,])

  }


  if(length(label)==numSeqs){
    featureMatrix<-as.data.frame(featureMatrix)
    featureMatrix<-cbind(featureMatrix,label)
  }

  row.names(featureMatrix)<-names(seqs)

  return(featureMatrix)
}



