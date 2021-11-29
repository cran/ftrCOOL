#' Z_curve_9bit_RNA (Zcurve9bit_RNA)
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
#' the number of columns is 9.
#'
#'
#' @export
#'
#' @examples
#'
#' fileLNC<-system.file("extdata/Carica_papaya101RNA.txt",package="ftrCOOL")
#' mat<-Zcurve9bit_RNA(seqs=fileLNC,ORF=TRUE,reverseORF=FALSE)


Zcurve9bit_RNA<-function(seqs,ORF=FALSE,reverseORF=TRUE,label=c()){


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
  featureMatrix<-matrix(0,nrow = numSeqs,ncol = 9)

  for(n in 1:numSeqs){
    seq<-seqs[n]
    chars<-unlist(strsplit(seq,NULL))
    lenSeq<-length(chars)
    PosInd1<-seq(1,lenSeq,3)
    PosInd2<-seq(2,lenSeq,3)
    PosInd3<-seq(3,lenSeq,3)

    charsPos1<-chars[PosInd1]
    charsPos2<-chars[PosInd2]
    charsPos3<-chars[PosInd3]

    tabCharsPos1<-table(charsPos1)
    tabCharsPos2<-table(charsPos2)
    tabCharsPos3<-table(charsPos3)


    pos1<-rep(0,4)
    names(pos1)<-c("A","C","G","T")
    pos1[names(tabCharsPos1)]<-tabCharsPos1


    pos2<-rep(0,4)
    names(pos2)<-c("A","C","G","U")
    pos2[names(tabCharsPos2)]<-tabCharsPos2


    pos3<-rep(0,4)
    names(pos3)<-c("A","C","G","U")
    pos3[names(tabCharsPos3)]<-tabCharsPos3



    x1=(pos1["A"]+pos1["G"])-(pos1["C"]+pos1["U"])
    y1=(pos1["A"]+pos1["C"])-(pos1["G"]+pos1["U"])
    z1=(pos1["A"]+pos1["U"])-(pos1["G"]+pos1["C"])

    x2=(pos2["A"]+pos2["G"])-(pos2["C"]+pos2["U"])
    y2=(pos2["A"]+pos2["C"])-(pos2["G"]+pos2["U"])
    z2=(pos2["A"]+pos2["U"])-(pos2["G"]+pos2["C"])

    x3=(pos3["A"]+pos3["G"])-(pos3["C"]+pos3["U"])
    y3=(pos3["A"]+pos3["C"])-(pos3["G"]+pos3["U"])
    z3=(pos3["A"]+pos3["U"])-(pos3["G"]+pos3["C"])

    vect<-c(x1,y1,z1,x2,y2,z2,x3,y3,z3)/lenSeq
    featureMatrix[n,]<-vect
  }
  colnames(featureMatrix)<-c("Pos_1.x","Pos_1.y","Pos_1.z","Pos_2.x","Pos_2.y","Pos_2.z","Pos_3.x","Pos_3.y","Pos_3.z")

  if(length(label)==numSeqs){
    featureMatrix<-as.data.frame(featureMatrix)
    featureMatrix<-cbind(featureMatrix,label)
  }

  return(featureMatrix)

}
