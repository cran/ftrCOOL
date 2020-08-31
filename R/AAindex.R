#' Amino Acid Index
#'
#' This function converts the amino acids of a sequence to a list of physicochemical properties in the aaIndex file.
#' For each amino acid, the function uses a numeric vector which shows the aaIndex of the amino acid.
#'
#'
#' @details In this function each amino acid is converted to a numeric vector. Elements of the vector represent a
#' physicochemical property for the amino acid.
#' In the aaIndex database, there are 554 amino acid indices. Users can choose the desired aaindex by specifying aaindexes through their ids or indexes in the aaIndex file, via selectedAAidx parameter.
#'
#' @note This function is provided for sequences with the same lengths.
#' Users can use 'txt' option in outFormat parameter for sequences with different lengths.
#' Warning: If outFormat is set to 'mat' for sequences with different lengths, it returns an error.
#' Also, when output format is 'txt', label information is not shown in the text file.
#' It is noteworthy that 'txt' format is not usable for machine learning purposes.
#'
#'
#' @param seqs is a FASTA file with amino acid sequences. Each sequence starts
#' with a '>' character. Also, seqs could be a string vector. Each element of the vector is a peptide/protein sequence.
#'
#'
#' @param selectedAAidx AAindex function works based on physicochemical properties. Users select the properties by their ids
#' or indexes in aaIndex2 file.
#'
#' @param standardized is a logical parameter. If it is set to TRUE, amino acid indices will be in the standard format.
#' The default value is TRUE.
#'
#'
#' @param threshold is a number between (0 , 1]. In selectedAAidx, indices with a correlation
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
#' @return The output depends on the outFormat parameter which can be either 'mat' or 'txt'. If outFormat is 'mat', the function returns a feature
#' matrix for sequences with the same length such that the number of columns is (sequence length)*(number of selected amino acid indexes)
#' and the number of rows is equal to the number of sequences. It is usable for machine learning purposes.
#' If the outFormat is 'txt', the output is written to a tab-delimited file.
#'
#'
#' @export
#'
#' @examples
#'
#' dir = tempdir()
#' ptmSeqsADR<-system.file("extdata/",package="ftrCOOL")
#' ptmSeqsVect<-as.vector(read.csv(paste0(ptmSeqsADR,"/ptmVect101AA.csv"))[,2])
#' mat<-AAindex(seqs = ptmSeqsVect, selectedAAidx=1:5,outFormat="mat")
#'
#' ad<-paste0(dir,"/aaidx.txt")
#' filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
#' AAindex(seqs = filePrs, selectedAAidx=1:5,standardized=TRUE,threshold=1,outFormat="txt"
#' ,outputFileDist=ad)
#'
#'

AAindex <- function(seqs,selectedAAidx=1:554,standardized=TRUE,threshold=1,label=c(),outFormat="mat",outputFileDist="")
{

  path.pack=system.file("extdata",package="ftrCOOL")


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
  lenSeqs<-sapply(seqs, nchar)

  numSeqs<-length(seqs)





  aaIdxAd<-paste0(path.pack,"/aaIndexcleaned_555.csv")
  aaIdx<- read.csv(aaIdxAd)
  if(standardized==TRUE){
    aaIdxAd<-paste0(path.pack,"/standardizedAAindex_555.csv")
    read.csv(aaIdxAd)

  }



  row.names(aaIdx)<-aaIdx[,1]
  aaIdx<-aaIdx[selectedAAidx,-1]
  aaIdx<-as.matrix(aaIdx)
  aaIdx<-type.convert(aaIdx)
  aaIdx<-aaIdx[-555,]

  ###deleteing high correlate features
  if(threshold<1){
    aaIdx<-t(aaIdx)
    corr<-cor(aaIdx)
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
      warMessage<-paste("The properties (",delPropName,") were deleted. They had a correlation with other properties more than the threshold")
      message(warMessage)
    }

    aaIdx<- aaIdx[,RMIDX]
    aaIdx<- t(aaIdx)
  }

  if(outFormat=="mat"){
    if(length(unique(lenSeqs))>1){
      stop("ERROR: All sequences should have the same length in 'mat' mode. For sequences with different lengths, please use 'txt' for outFormat parameter")
    }
    dimFeaVct <-nrow(aaIdx)*lenSeqs[1]
    featureMatrix=matrix(0,nrow = numSeqs,ncol = dimFeaVct)


    tempname<-rep(1:lenSeqs[1],each=nrow(aaIdx))
    temp2name<-rep(row.names(aaIdx),lenSeqs[1])
    colna<-paste0(temp2name,"-pos",tempname)
    colnames(featureMatrix)<-colna

    for(n in 1:numSeqs){
      seq<-seqs[n]
      chars<-unlist(strsplit(seq,NULL))
      aaMatselected<-aaIdx[,chars]
      vect<-as.vector(aaMatselected)
      featureMatrix[n,]<-t(vect)


    }
    if(length(label)==numSeqs){
      featureMatrix<-as.data.frame(featureMatrix)
      featureMatrix<-cbind(featureMatrix,label)
    }
    row.names(featureMatrix)<-names(seqs)
    return(featureMatrix)
  }
  else{
    nameSeq<-names(seqs)
    for(n in 1:numSeqs){
      seq<-seqs[n]
      chars<-unlist(strsplit(seq,split = ""))
      aaMatselected<-aaIdx[,chars]
      aaMatselected<-as.character(aaMatselected)
      temp<-c(nameSeq[n],aaMatselected)
      temp<-paste(temp,collapse = "\t")
      write(temp,outputFileDist,append = TRUE)

    }
  }
}
