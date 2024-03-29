#' Trinucleotide physicochemical properties (TPCP_DNA)
#'
#' This function replaces trinucleotides in a sequence with their physicochemical properties which is multiplied by normalized frequency of that tri-nucleotide.
#'
#'
#' @details There are 12 physicochemical indexes in the trinucleotide database.
#'
#' @note This function is provided for sequences with the same lengths.
#' Users can use 'txt' option in outFormat for sequences with different lengths.
#' Warning: If outFormat is set to 'mat' for sequences with different lengths, it returns an error.
#' Also, when output format is 'txt', label information is not shown in the text file.
#' It is noteworthy that 'txt' format is not usable for machine learning purposes if sequences have different sizes. Otherwise 'txt' format
#' is also usable for machine learning purposes.
#'
#' @param seqs is a FASTA file containing nucleotide sequences. The sequences start
#' with '>'. Also, seqs could be a string vector. Each element of the vector is a nucleotide sequence.
#'
#' @param selectedIdx TPCP_DNA function works based on physicochemical properties. Users, select the properties by their ids
#' or indexes in TRI_DNA index file.
#' The default values of the vector are the  ids in "Dnase I", "Bendability (DNAse)".
#'
#'
#' @param threshold is a number between 0 to 1. In selectedIdx, indices with a correlation
#' higher than the threshold will be deleted. The default value is 1.
#'
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @param outFormat (output format) can take two values: 'mat'(matrix) and 'txt'. The default value is 'mat'.
#'
#' @param outputFileDist shows the path and name of the 'txt' output file.
#'
#'
#' @return The output depends on the outFormat parameter which can be either 'mat' or 'txt'. If outFormat is 'mat', the function returns a feature
#' matrix for sequences with the same length such that the number of columns is (sequence length-2)*(number of selected trinucleotide properties)
#' and the number of rows is equal to the number of sequences.
#' If the outFormat is 'txt', the output is written to a tab-delimited file.
#'
#'
#' @export
#'
#' @examples
#'
#' fileLNC<-system.file("extdata/Athaliana1.fa",package="ftrCOOL")
#' vect<-TPCP_DNA(seqs = fileLNC,threshold=1,outFormat="mat")



TPCP_DNA <- function(seqs,selectedIdx=c("Dnase I", "Bendability (DNAse)"),threshold=1,label=c(),outFormat="mat",outputFileDist="")
{

  path.pack=system.file("extdata",package="ftrCOOL")
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

  lenSeqs<-sapply(seqs, nchar)

  numSeqs<-length(seqs)



  nucIdxAd<-paste0(path.pack,"/TRI_DNA.csv")
  #nucIdxAd<-paste0(path.pack,"/DI_DNA.csv")

  nucIdx<-read.csv(nucIdxAd)

  row.names(nucIdx)<-nucIdx[,1]
  nucIdx<-nucIdx[selectedIdx,-1]
  nucIdx<-as.matrix(nucIdx)
  nucIdx<-type.convert(nucIdx)




  ###deleteing high correlate features
  if(threshold<1){
    nucIdx<-t(nucIdx)
    corr<-cor(nucIdx)
    corr2<-corr^2

    tmp<-corr2
    tmp[upper.tri(tmp)]<-0

    for(i in 1:ncol(tmp)){
      tmp[i,i]=0
    }

    TFidx<-apply(tmp,2,function(x) any(x > threshold))
    DlIDX<-which(TFidx==TRUE)
    RMIDX<-which(TFidx==FALSE)


    if(length(DlIDX)!=0){
      delPropName<-names(DlIDX)
      delPropName<-toString(delPropName)
      warMessage<-paste("The properties (",delPropName,") were deleted. They had a correlation with other properties more than the threshold ")
      message(warMessage)
    }

    nucIdx<- nucIdx[,RMIDX]
    nucIdx<- t(nucIdx)
  }



  dimFeaVct <-nrow(nucIdx)*64
  featureMatrix=matrix(0,nrow = numSeqs,ncol = dimFeaVct)
  temp1<-rep(rownames(nucIdx),each=64)
  temp2<-rep(colnames(nucIdx),nrow(nucIdx))
  colnames(featureMatrix)<-paste0(temp1,"_",temp2)


  for(n in 1:numSeqs){
    seq<-seqs[n]
    chars<-unlist(strsplit(seq,NULL))
    lenSeq<-length(chars)
    temp1<-chars[1:(lenSeq-2)]
    temp2<-chars[2:(lenSeq-1)]
    temp3<-chars[3:lenSeq]
    triMer<-paste0(temp1,temp2,temp3,"")
    freqTriMer<-table(triMer)


    a=apply(nucIdx,MARGIN = 1,function(x){x[names(freqTriMer)]*freqTriMer/(lenSeq-2)
    })

    vect<-as.vector(a)
    tempN1<-rep(colnames(a),each=nrow(a))
    tempN2<-rep(rownames(a),ncol(a))
    names(vect)<-paste0(tempN1,"_",tempN2)
    featureMatrix[n,names(vect)]<-vect


  }

  if(length(label)==numSeqs){
    featureMatrix<-as.data.frame(featureMatrix)
    featureMatrix<-cbind(featureMatrix,label)
  }
  row.names(featureMatrix)<-names(seqs)
  return(featureMatrix)


}
