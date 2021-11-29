#' Z_curve_144bit_RNA (Zcurve144bit_RNA)
#'
#' These group of functions (Zcurve (9, 12, 36, 48, 144)_bit) function calculates the Z-curves. Z-curves are based on freqiencies of ribonucleotides, di-ribonucleotides, or tri-ribonucleotides and their positions on the sequences.
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
#' the number of columns is 144.
#'
#' @export
#'
#' @examples
#'
#' fileLNC<-system.file("extdata/Carica_papaya101RNA.txt",package="ftrCOOL")
#' mat<-Zcurve144bit_RNA(seqs=fileLNC,ORF=TRUE,reverseORF=FALSE)


Zcurve144bit_RNA<-function(seqs,ORF=FALSE,reverseORF=TRUE,label=c()){


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

  er<-which(lenSeqs<5)
  if(length(er)!=0){
    stop(paste("Sequnse/Sequnses",er,"is/are smaller than minSize(5)"))
  }


  DiNuc<-nameKmer(k=2,type = "rna")
  TriNuc<-nameKmer(k=3,type = "rna")
  featureMatrix<-matrix(ncol = 144,nrow = numSeqs)
  for(n in 1:numSeqs){

    seq<-seqs[n]
    lenSeq<-lenSeqs[n]



    PosInd1<-seq(1,lenSeq,3)
    PosInd1<-PosInd1[(PosInd1+2)<=lenSeq]
    EndPos1<-PosInd1+2

    PosInd2<-seq(2,lenSeq,3)
    PosInd2<-PosInd2[(PosInd2+2)<=lenSeq]
    EndPos2<-PosInd2+2


    PosInd3<-seq(3,lenSeq,3)
    PosInd3<-PosInd3[(PosInd3+2)<=lenSeq]
    EndPos3<-PosInd3+2

    TrimPos1<-substring(seq,PosInd1,EndPos1)
    TrimPos2<-substring(seq,PosInd2,EndPos2)
    TrimPos3<-substring(seq,PosInd3,EndPos3)

    freqTrimPos1<-rep(0,64)
    names(freqTrimPos1)<-TriNuc

    freqTrimPos2<-rep(0,64)
    names(freqTrimPos2)<-TriNuc

    freqTrimPos3<-rep(0,64)
    names(freqTrimPos3)<-TriNuc

    TrimPos1<-table(TrimPos1)
    TrimPos2<-table(TrimPos2)
    TrimPos3<-table(TrimPos3)

    freqTrimPos1[names(TrimPos1)]<-TrimPos1
    freqTrimPos2[names(TrimPos2)]<-TrimPos2
    freqTrimPos3[names(TrimPos3)]<-TrimPos3



    Pxa=freqTrimPos1[paste(DiNuc,"A",sep="")]
    Pxg=freqTrimPos1[paste(DiNuc,"G",sep="")]
    Pxc=freqTrimPos1[paste(DiNuc,"C",sep="")]
    Pxt=freqTrimPos1[paste(DiNuc,"U",sep="")]
    xVect1<-(Pxa+Pxg)-(Pxc+Pxt)
    yVect1<-(Pxa+Pxc)-(Pxg+Pxt)
    zVect1<-(Pxa+Pxt)-(Pxg+Pxc)


    Pxa2=freqTrimPos2[paste(DiNuc,"A",sep="")]
    Pxg2=freqTrimPos2[paste(DiNuc,"G",sep="")]
    Pxc2=freqTrimPos2[paste(DiNuc,"C",sep="")]
    Pxt2=freqTrimPos2[paste(DiNuc,"U",sep="")]
    xVect2<-(Pxa2+Pxg2)-(Pxc2+Pxt2)
    yVect2<-(Pxa2+Pxc2)-(Pxg2+Pxt2)
    zVect2<-(Pxa2+Pxt2)-(Pxg2+Pxc2)




    Pxa3=freqTrimPos3[paste(DiNuc,"A",sep="")]
    Pxg3=freqTrimPos3[paste(DiNuc,"G",sep="")]
    Pxc3=freqTrimPos3[paste(DiNuc,"C",sep="")]
    Pxt3=freqTrimPos3[paste(DiNuc,"U",sep="")]
    xVect3<-(Pxa3+Pxg3)-(Pxc3+Pxt3)
    yVect3<-(Pxa3+Pxc3)-(Pxg3+Pxt3)
    zVect3<-(Pxa3+Pxt3)-(Pxg+Pxc3)



    mat<-rbind(xVect1,yVect1,zVect1,xVect2,yVect2,zVect2,xVect3,yVect3,zVect3)
    featureMatrix[n,]<-as.vector(mat)/(lenSeq-2)

  }

  tempName1<-rep(DiNuc,each=9)
  tempName2<-rep(c("x","y","z"),48)
  tempName3<-rep(c("1","2","3"),each=3)
  tmp<-paste0("Pos.",tempName3,".",tempName1,".",tempName2)
  colnames(featureMatrix)<-tmp


  if(length(label)==numSeqs){
    featureMatrix<-as.data.frame(featureMatrix)
    featureMatrix<-cbind(featureMatrix,label)
  }

  return(featureMatrix)
}
