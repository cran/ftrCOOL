#' Enhanced nucleotide Composition
#'
#' This function slides a window over the input sequence(s).
#' Also, it computes the composition of nucleotides that appears within the limits of the window.
#' Please note that when a feature vector of sequences is given as the input, their length should be equal.
#' Otherwise, either an error (in case the outFortmat is 'mat') or a text file (when the outFortmat is 'txt') is returned.
#' Please note that the text file is not suitable for machine learning purposes.
#'
#' @note When overlap is FALSE, the last partition represented by the window may have a different length with other parts.
#'
#' @param seqs is a FASTA file containing nucleotide sequences. The sequences start
#' with '>'. Also, seqs could be a string vector. Each element of the vector is a nucleotide sequence.
#'
#' @param winSize is a number which shows the size of the window.
#'
#' @param overLap This parameter shows how the window
#' moves on the sequence. If the overlap is set to TRUE, the next window would have distance 1 with
#' the previous window. Otherwise, the next window will start from the next amino acid after the previous window.
#' There is no overlap between the next and previous windows.
#'
#'
#' @param outFormat (output format) can take two values: 'mat'(matrix) and 'txt'. The default value is 'mat'.
#'
#' @param outputFileDist shows the path and name of the 'txt' output file.
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @return The output depends on the outFormat parameter which can be either 'mat' or 'txt'. If outFormat is 'mat', the function returns a feature
#' matrix for sequences with the same length such that the number of columns is (4 * number of partitions displayed by the window)
#' and the number of rows is equal to the number of sequences. It is usable for machine learning purposes.
#' If the outFormat is 'txt', the output is written to a tab-delimited file.
#'
#' @export
#'
#' @examples
#'
#' dir = tempdir()
#' LNCSeqsADR<-system.file("extdata/",package="ftrCOOL")
#' LNC50Nuc<-as.vector(read.csv(paste0(LNCSeqsADR,"/LNC50Nuc.csv"))[,2])
#' mat<-ENUComposition(seqs = LNC50Nuc, winSize=20,outFormat="mat")
#'
#' ad<-paste0(dir,"/ENUCcompos.txt")
#' fileLNC<-system.file("extdata/Athaliana_LNCRNA.fa",package="ftrCOOL")
#' ENUComposition(seqs = fileLNC,outFormat="txt",winSize=20
#' ,outputFileDist=ad,overLap=FALSE)


ENUComposition <- function(seqs,winSize=20,overLap=TRUE,outFormat='mat',outputFileDist="",label=c())
{

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

  lenSeqs<-sapply(seqs,nchar)

  if(!all(lenSeqs>=winSize)){
    deletInd<-which(lenSeqs<winSize)
    deletedNames<-names(deletInd)
    strNames<-toString(deletedNames)
    lens<-lenSeqs[deletInd]
    strlens<-toString(lens)

    warning(paste("Sequences",strNames,"with lengths",strlens,"were deleted. Their lenghts were smaller than the window size"))

    if(length(label)==length(lenSeqs)){
      label<-label[-deletInd]
    }
    lenSeqs<-lenSeqs[-deletInd]
    seqs<-seqs[-deletInd]
  }
  dict<-list("A"=1,"C"=2,"G"=3,"T"=4)


  numSeqs<-length(seqs)
  flag=1
  if(overLap==FALSE){
    st<-0:(ceiling(lenSeqs[1]/winSize)-1)
    st<-(st*winSize)+1
    en<-st+winSize-1
    len=ceiling(lenSeqs[1]/winSize)
    flag=0
  }
  else{
    st<-seq(1,(lenSeqs[1]-winSize+1),1)
    en<-seq(winSize,lenSeqs[1],1)
    len<-lenSeqs[1]-winSize+1
  }

  if(outFormat=='mat'){

  if(length(unique(lenSeqs))>1){
    stop("ERROR: All sequences should have the same length in 'mat' mode. For sequences with different lengths, please use 'txt' for outFormat parameter")
  }




  featureMatrix<-matrix(0,ncol = (4*len),nrow = numSeqs)

  tname<-nameKmer(k=1,type = "dna")
  tname<-rep(tname,len)
  wname<-rep(1:len,each =4)

  coln<-paste(tname,"w",wname,sep = "")

  colnames(featureMatrix)<-coln


  for(n in 1:numSeqs){

    seq<-seqs[n]
    subSeqs<-substring(seq,st,en)
    chars<-lapply(subSeqs, function(i) unlist(strsplit(i,"")))


    k<-length(chars)
    for(j in 1:k)
    {
      twinSeq<-table(chars[[j]])
      ntwinseq<-names(twinSeq)

      for(g in 1:length(twinSeq))
      {
        colindex<-(j-1)*4+as.numeric(dict[ntwinseq[g]])
        featureMatrix[n,colindex]<-twinSeq[g]

      }
    }


  }


  if(length(label)==numSeqs){
    featureMatrix<-as.data.frame(featureMatrix)
    featureMatrix<-cbind(featureMatrix,label)
  }
  row.names(featureMatrix)<-names(seqs)
  return(featureMatrix)


  }



  else if(outFormat=="txt"){

    nameSeq<-names(seqs)

    for(n in 1:numSeqs){

      seq<-seqs[n]
      subSeqs<-substring(seq,st,en)
      chars<-lapply(subSeqs, function(i) unlist(strsplit(i,"")))


      k<-length(chars)
      vect<-vector(mode = "numeric",length = 4*len)
      for(j in 1:k)
      {
        twinSeq<-table(chars[[j]])
        ntwinseq<-names(twinSeq)

        for(g in 1:length(twinSeq))
        {
          colindex<-(j-1)*4+as.numeric(dict[ntwinseq[g]])
          vect[colindex]<-twinSeq[g]

        }
      }

      temp<-c(nameSeq[n],vect)
      temp2<-paste(temp,collapse = "\t")
      write(temp2, outputFileDist, append=TRUE)

    }
    }

}



