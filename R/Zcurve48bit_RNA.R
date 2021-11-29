#' Z_curve_48bit_RNA (Zcurve48bit_RNA)
#'
#' These group of functions (Zcurve (9, 12, 36, 48, 144)_bit) function calculates the Z-curves. Z-curves are based on freqiencies of ribo
#' ribonucleotides, di-ribonucleotides, or tri-ribonucleotides and their positions on the sequences.
#' For more information about the methods please refer to reference part.
#'
#' @references Gao,F. and Zhang,C.T. Comparison of various algorithms for recognizing short coding sequences of human genes. Bioinformatics, (2004).
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
#' @return This function returns a feature matrix. The number of rows is equal to the number of sequences and
#' the number of columns is 48.
#'
#'
#' @export
#'
#' @examples
#'
#' fileLNC<-system.file("extdata/Carica_papaya101RNA.txt",package="ftrCOOL")
#' mat<-Zcurve48bit_RNA(seqs=fileLNC,ORF=TRUE,reverseORF=FALSE)


Zcurve48bit_RNA<-function(seqs,ORF=FALSE,reverseORF=TRUE,label=c()){


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

  TriNuc<-nameKmer(k=3,type = "rna")
  featureMatrix<-matrix(ncol = 48,nrow = numSeqs)
  for(n in 1:numSeqs){
    seq<-seqs[n]
    lenSeq<-lenSeqs[n]
    seqChars<-unlist(strsplit(seq,split = ""))
    temp<-seqChars[1:(lenSeq-2)]
    temp2<-seqChars[2:(lenSeq-1)]
    temp3<-seqChars[3:lenSeq]
    Trimer<-paste(temp,temp2,temp3,sep = "")
    tabtrimer<-table(Trimer)
    freqTrim<-rep(0,64)
    names(freqTrim)<-TriNuc
    freqTrim[names(tabtrimer)]<-tabtrimer

    diNucs<-nameKmer(k=2,type = "rna")


    Pxa=freqTrim[paste(diNucs,"A",sep="")]
    Pxg=freqTrim[paste(diNucs,"G",sep="")]
    Pxc=freqTrim[paste(diNucs,"C",sep="")]
    Pxt=freqTrim[paste(diNucs,"U",sep="")]
    xVect<-(Pxa+Pxg)-(Pxc+Pxt)
    yVect<-(Pxa+Pxc)-(Pxg+Pxt)
    zVect<-(Pxa+Pxt)-(Pxg+Pxc)

    mat<-rbind(xVect,yVect,zVect)

    vecMat<-as.vector(mat)

    featureMatrix[n,]<-vecMat/(lenSeq-2)

  }
  tempName1<-rep(diNucs,each=3)
  tempName2<-rep(c("x","y","z"),16)
  tmp<-paste0(tempName1,".",tempName2)
  colnames(featureMatrix)<-tmp


  if(length(label)==numSeqs){
    featureMatrix<-as.data.frame(featureMatrix)
    featureMatrix<-cbind(featureMatrix,label)
  }

  return(featureMatrix)

}
