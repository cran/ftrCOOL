#' Fickett Score (fickettScore)
#'
#' This function calculates the ficket score of each sequence.
#'
#' @param seqs is a FASTA file containing nucleotide sequences. The sequences start
#' with '>'. Also, seqs could be a string vector. Each element of the vector is a nucleotide sequence.
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
#' @return The function returns a feature vector. The length of the vector is equal to the number of sequences.
#' Each entry in the vector contains the value of the fickett score.
#'
#'
#' @export
#'
#' @examples
#'
#' fileLNC<-system.file("extdata/Athaliana_LNCRNA.fa",package="ftrCOOL")
#' vect<-fickettScore(seqs=fileLNC,ORF=TRUE,reverseORF=FALSE)


fickettScore<-function(seqs,ORF=FALSE,reverseORF=TRUE,label=c()){

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

  numSeqs<-length(seqs)
  lenSeqs<-sapply(seqs,nchar)

  fickettScore<-vector(mode = "numeric",length = numSeqs)
  positionProb = list("A" = c(0.94, 0.68, 0.84, 0.93, 0.58, 0.68,
                              0.45, 0.34, 0.2, 0.22), "C" = c(0.8, 0.7, 0.7, 0.81, 0.66,
                                                              0.48, 0.51, 0.33, 0.3, 0.23), "G" = c(0.9, 0.88, 0.74,
                                                                                                    0.64, 0.53, 0.48, 0.27, 0.16, 0.08, 0.08), "T" = c(0.97,
                                                                                                                                                       0.97, 0.91, 0.68, 0.69, 0.44, 0.54, 0.2, 0.09, 0.09))
  positionWeight = c("A" = 0.26, "C" = 0.18, "G" = 0.31, "T" = 0.33)
  positionParam = c(1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 0)

  contentProb = list("A" = c(0.28, 0.49, 0.44, 0.55, 0.62, 0.49,
                             0.67, 0.65, 0.81, 0.21), "C" = c(0.82, 0.64, 0.51, 0.64,
                                                              0.59, 0.59, 0.43, 0.44, 0.39, 0.31), "G" = c(0.4, 0.54,
                                                                                                           0.47, 0.64, 0.64, 0.73, 0.41, 0.41, 0.33, 0.29), "T" = c(0.28,
                                                                                                                                                                    0.24, 0.39, 0.4, 0.55, 0.75, 0.56, 0.69, 0.51, 0.58))
  contentWeight = c("A" = 0.11, "C" = 0.12, "G" = 0.15, "T" = 0.14)
  contentParam = c(0.33, 0.31, 0.29, 0.27, 0.25, 0.23, 0.21,0.19, 0.17, 0)


  for(n in 1:numSeqs){
    seq<-seqs[n]
    chars<-unlist(strsplit(seqs[n],NULL))
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

    pos1<-tabCharsPos1[c("A","C","G","T")]
    pos1[is.na(pos1)]<-0

    pos2<-tabCharsPos2[c("A","C","G","T")]
    pos2[is.na(pos2)]<-0

    pos3<-tabCharsPos3[c("A","C","G","T")]
    pos3[is.na(pos3)]<-0

    Apos<-max(pos1[1],pos2[1],pos3[1])/(min(pos1[1],pos2[1],pos3[1])+1)
    Cpos<-max(pos1[2],pos2[2],pos3[2])/(min(pos1[2],pos2[2],pos3[2])+1)
    Gpos<-max(pos1[3],pos2[3],pos3[3])/(min(pos1[3],pos2[3],pos3[3])+1)
    Tpos<-max(pos1[4],pos2[4],pos3[4])/(min(pos1[4],pos2[4],pos3[4])+1)

    tabWholeSeq<-table(chars)
    tabWholeSeq<-tabWholeSeq[c("A","C","G","T")]
    tabWholeSeq[is.na(tabWholeSeq)]<-0

    Acontent<-tabWholeSeq[1]/lenSeq
    Ccontent<-tabWholeSeq[2]/lenSeq
    Gcontent<-tabWholeSeq[3]/lenSeq
    Tcontent<-tabWholeSeq[4]/lenSeq



    indA<-which(positionParam<=Apos)
    AprobPos<-positionProb[["A"]][indA[1]]

    indC<-which(positionParam<=Cpos)
    CprobPos<-positionProb[["C"]][indC[1]]



    indG<-which(positionParam<=Gpos)
    GprobPos<-positionProb[["G"]][indG[1]]

    indT<-which(positionParam<=Tpos)
    TprobPos<-positionProb[["T"]][indT[1]]


    indA<-which(contentParam<=Acontent)
    AprobCon<-contentProb[["A"]][indA[1]]


    indC<-which(contentParam<=Ccontent)
    CprobCon<-contentProb[["C"]][indC[1]]

    indG<-which(contentParam<=Gcontent)
    GprobCon<-contentProb[["G"]][indG[1]]

    indT<-which(contentParam<=Tcontent)
    TprobCon<-contentProb[["T"]][indT[1]]

    fickettScore[n]<-((AprobPos*positionWeight["A"])+(CprobPos*positionWeight["C"])+
                        (GprobPos*positionWeight["G"])+(TprobPos*positionWeight["T"])+
                        (AprobCon*contentWeight["A"])+(CprobCon*contentWeight["C"])+
                        (GprobCon*contentWeight["G"])+(TprobCon*contentWeight["T"]))



  }

  names(fickettScore)<-names(seqs)
  if(length(label)==numSeqs){

    fickettScore<-cbind(fickettScore,label)
  }

  return(fickettScore)

}
