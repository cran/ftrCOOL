#' Codon Usage in RNA (CodonUsage_RNA)
#'
#' This function calculates the codon usage for each sequence.
#'
#' @param seqs is a FASTA file containing ribonucleotide sequences. The sequences start
#' with '>'. Also, seqs could be a string vector. Each element of the vector is a ribonucleotide sequence.
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
#'
#' @return A feature matrix such that the number of columns is 4^3 and the number of rows is equal to the number of sequences.
#'
#'
#' @export
#'
#' @examples
#'
#' fileLNC<-system.file("extdata/Carica_papaya101RNA.txt",package="ftrCOOL")
#' mat<-CodonUsage_RNA(seqs=fileLNC,ORF=TRUE,reverseORF=FALSE)
#'


CodonUsage_RNA<-function(seqs,ORF=FALSE,reverseORF=TRUE,label=c())
{

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

  dict<-list("A"=1,"C"=2,"G"=3,"U"=4)


  flag=0
  if(ORF==TRUE){
    if(length(label)==length(seqs)){
      names(label)=names(seqs)
      flag=1
    }
    seqs=maxORF_RNA(seqs,reverse=reverseORF)
    if(flag==1)
      label=label[names(seqs)]
  }
  numSeqs<-length(seqs)
  lenSeqs<-sapply(seqs,nchar)
  featureMatrix<-matrix(0,nrow = numSeqs,ncol = 4^3)
  colnames(featureMatrix)<-paste0(nameKmer(k=3,type = "rna"),"(usg)")
  l=3

  for(n in 1:numSeqs){
    seq<-seqs[n]
    triplet<-substring(seq,seq(1,(lenSeqs[n]-2),3),seq(3,lenSeqs[n],3))
    triplet<-triplet[1:(floor(lenSeqs[n]/3))]
    tableContent<-table(triplet)
    tabNames<-names(tableContent)

    for(i in 1:length(tableContent))
    {
      temp<-unlist(strsplit(tabNames[i],split = ""))
      num=0
      for(j in 1:l){
        pow<-4^(l-j)
        num<-num+(((as.numeric(dict[temp[j]]))-1)*pow)
      }
      num<-num+1
      featureMatrix[n,num]<-tableContent[i]
    }

  }
  if(length(label)==numSeqs){
    featureMatrix<-as.data.frame(featureMatrix)
    featureMatrix<-cbind(featureMatrix,label)
  }
  row.names(featureMatrix)<-names(seqs)
  return(featureMatrix)

}
