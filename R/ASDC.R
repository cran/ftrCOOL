#' Adaptive skip dipeptide composition (ASDC)
#'
#' This descriptor sufficiently considers the correlation information present not only between adjacent residues but also between intervening residues.
#' This function calculates frequency of pair amino acids omitting gaps between them. Then this function normalizes each value through dividing each frequency by summition(frequencies).
#'
#'@references Wei L, Zhou C, Chen H, Song J, Su R. ACPred-FL: a sequence-based predictor using effective feature representation to improve the prediction of anti-cancer peptides. Bioinformatics (2018).
#'
#' @param seqs is a FASTA file with amino acid sequences. Each sequence starts
#' with a '>' character. Also, seqs could be a string vector. Each element of the vector is a peptide/protein sequence.
#'
#'
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @return The function returns a feature matrix. The number of rows is equal to the number of sequences and
#' the number of columns is 400 (all posible amino acid pairs).
#'
#'
#' @export
#'
#' @examples
#'
#'
#' filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
#' mat<-ASDC(seqs=filePrs)
#'



ASDC<- function(seqs,label=c()){


  upto=TRUE


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
  flag=0


  numSeqs=length(seqs)

  rng<-sapply(seqs, nchar)
  rng<-rng-1





  featureMatrix <- matrix(0 , ncol = 400,nrow = numSeqs)
  dipep<-nameKmer(k=2,type = "aa")
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


      tempvect<-vector(mode = "numeric",length = 400)
      names(tempvect)<-tempname
      tempvect[nmtbkmers]<-tbkmers
      featureMatrix[n,]<-featureMatrix[n,]+tempvect
    }
    featureMatrix[n,]<-featureMatrix[n,]/sum(featureMatrix[n,])

  }

  # if(normalized==TRUE){
  #   seqLen<-sapply(seqs, nchar)
  #   #featureMatrix<-featureMatrix/seqLen
  # }
  if(length(label)==numSeqs){
    featureMatrix<-as.data.frame(featureMatrix)
    featureMatrix<-cbind(featureMatrix,label)
  }

  row.names(featureMatrix)<-names(seqs)

  return(featureMatrix)
}



