#' Z_curve_36bit_RNA (Zcurve36bit_RNA)
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
#' the number of columns is 36.
#'
#' @export
#'
#' @examples
#'
#' fileLNC<-system.file("extdata/Carica_papaya101RNA.txt",package="ftrCOOL")
#' mat<-Zcurve36bit_RNA(seqs=fileLNC,ORF=TRUE,reverseORF=FALSE)


Zcurve36bit_RNA<-function(seqs,ORF=FALSE,reverseORF=TRUE,label=c()){


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



  DiNuc<-nameKmer(k=2,type = "rna")
  featureMatrix<-matrix(ncol = 36,nrow = numSeqs)
  for(n in 1:numSeqs){

    seq<-seqs[n]
    lenSeq<-lenSeqs[n]



    PosInd1<-seq(1,lenSeq,3)
    PosInd1<-PosInd1[(PosInd1+1)<=lenSeq]
    EndPos1<-PosInd1+1

    PosInd2<-seq(2,lenSeq,3)
    PosInd2<-PosInd2[(PosInd2+1)<=lenSeq]
    EndPos2<-PosInd2+1


    PosInd3<-seq(3,lenSeq,3)
    PosInd3<-PosInd3[(PosInd3+1)<=lenSeq]
    EndPos3<-PosInd3+1


    DimPos1<-substring(seq,PosInd1,EndPos1)
    DimPos2<-substring(seq,PosInd2,EndPos2)
    DimPos3<-substring(seq,PosInd3,EndPos3)

    freqDimPos1<-rep(0,16)
    names(freqDimPos1)<-DiNuc

    freqDimPos2<-rep(0,16)
    names(freqDimPos2)<-DiNuc

    freqDimPos3<-rep(0,16)
    names(freqDimPos3)<-DiNuc

    DimPos1<-table(DimPos1)
    DimPos2<-table(DimPos2)
    DimPos3<-table(DimPos3)

    freqDimPos1[names(DimPos1)]<-DimPos1
    freqDimPos2[names(DimPos2)]<-DimPos2
    freqDimPos3[names(DimPos3)]<-DimPos3




    Nucs<-c("A","C","G","U")
    mat1<-matrix(nrow = 3,ncol = 4)
    rownames(mat1)<-c("x","y","z")
    colnames(mat1)<-c("A","C","G","U")
    for(i in 1:4){
      Pxa=freqDimPos1[paste(Nucs[i],"A",sep="")]
      Pxg=freqDimPos1[paste(Nucs[i],"G",sep="")]
      Pxc=freqDimPos1[paste(Nucs[i],"C",sep="")]
      Pxt=freqDimPos1[paste(Nucs[i],"U",sep="")]
      mat1["x",i]<-(Pxa+Pxg)-(Pxc+Pxt)
      mat1["y",i]<-(Pxa+Pxc)-(Pxg+Pxt)
      mat1["z",i]<-(Pxa+Pxt)-(Pxg+Pxc)
    }

    mat2<-matrix(nrow = 3,ncol = 4)
    rownames(mat2)<-c("x","y","z")
    colnames(mat2)<-c("A","C","G","U")
    for(i in 1:4){
      Pxa=freqDimPos2[paste(Nucs[i],"A",sep="")]
      Pxg=freqDimPos2[paste(Nucs[i],"G",sep="")]
      Pxc=freqDimPos2[paste(Nucs[i],"C",sep="")]
      Pxt=freqDimPos2[paste(Nucs[i],"U",sep="")]
      mat2["x",i]<-(Pxa+Pxg)-(Pxc+Pxt)
      mat2["y",i]<-(Pxa+Pxc)-(Pxg+Pxt)
      mat2["z",i]<-(Pxa+Pxt)-(Pxg+Pxc)
    }

    mat3<-matrix(nrow = 3,ncol = 4)
    rownames(mat3)<-c("x","y","z")
    colnames(mat3)<-c("A","C","G","U")
    for(i in 1:4){
      Pxa=freqDimPos3[paste(Nucs[i],"A",sep="")]
      Pxg=freqDimPos3[paste(Nucs[i],"G",sep="")]
      Pxc=freqDimPos3[paste(Nucs[i],"C",sep="")]
      Pxt=freqDimPos3[paste(Nucs[i],"U",sep="")]
      mat3["x",i]<-(Pxa+Pxg)-(Pxc+Pxt)
      mat3["y",i]<-(Pxa+Pxc)-(Pxg+Pxt)
      mat3["z",i]<-(Pxa+Pxt)-(Pxg+Pxc)
    }


    mat<-rbind(mat1,mat2,mat3)
    featureMatrix[n,]<-as.vector(mat)/(lenSeq-1)


  }

  tempName1<-rep(c("A","C","G","U"),each=9)
  tempName2<-rep(c("x","y","z"),12)
  temp3<-rep(c("1","2","3"),each=3)
  colnames(featureMatrix)<-paste0("Pos.",temp3,".",tempName1,".",tempName2)

  if(length(label)==numSeqs){
    featureMatrix<-as.data.frame(featureMatrix)
    featureMatrix<-cbind(featureMatrix,label)
  }

  return(featureMatrix)
}
