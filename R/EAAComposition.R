#' Enhanced Amino Acid Composition (EAAComposition)
#'
#' This function slides a window over the input sequence(s).
#' Also, it computes the composition of amino acids that appears within the limits of the window.
#'
#' @note This function is provided for sequences with the same lengths.
#' Users can use 'txt' option in outFormat for sequences with different lengths.
#' Warning: If outFormat is set to 'mat' for sequences with different lengths, it returns an error.
#' Also, when output format is 'txt', label information is not shown in the text file.
#' It is noteworthy that 'txt' format is not usable for machine learning purposes if sequences have different sizes. Otherwise 'txt' format
#' is also usable for machine learning purposes.
#'
#' @references Chen, Zhen, et al. "iFeature: a python package and web server for features extraction and selection from protein and peptide sequences." Bioinformatics 34.14 (2018): 2499-2502.
#'
#' @details Column names in the output matrix are Wi(aa), where aa shows an amino acid type ("A", "C", "D",..., "Y") and
#' i indicates the number of times that the window has moved over the sequence(s).
#'
#' @note When overlap is FALSE, the last partition represented by the window may have a different length with other parts.
#'
#'
#' @param seqs is a FASTA file with amino acid sequences. Each sequence starts
#' with a '>' character. Also, seqs could be a string vector. Each element of the vector is a peptide/protein sequence.
#'
#' @param winSize is a number which shows the size of the window.
#'
#' @param overLap This parameter shows how the window
#' moves over the sequence. If overlap is set to FALSE, the window slides over the sequence in such a way that every time the window moves, it covers a unique portion of the sequence.
#' Otherwise, portions of the sequence which appear within the window limits have "winSize-1" amino acids in common.
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
#' matrix for sequences with the same length such that the number of columns is (20 * number of partitions displayed by the window)
#' and the number of rows is equal to the number of sequences.
#' If the outFormat is 'txt', the output is written to a tab-delimited file.
#'
#'
#'
#' @export
#'
#' @examples
#'
#' dir = tempdir()
#' ptmSeqsADR<-system.file("extdata/",package="ftrCOOL")
#' ptmSeqsVect<-as.vector(read.csv(paste0(ptmSeqsADR,"/ptmVect101AA.csv"))[,2])
#' mat<-EAAComposition(seqs = ptmSeqsVect,winSize=50, overLap=FALSE,outFormat='mat')
#'
#' ad<-paste0(dir,"/EaaCompos.txt")
#' filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
#' EAAComposition(seqs = filePrs,winSize=50, overLap=FALSE,outFormat="txt"
#' ,outputFileDist=ad)
#'
#' unlink("dir", recursive = TRUE)





EAAComposition <- function(seqs,winSize=50,overLap=TRUE,label=c(),outFormat='mat',outputFileDist="")
{

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

  lenSeqs<-sapply(seqs,nchar)
  if(!all(lenSeqs>=winSize)){
    deletInd<-which(lenSeqs<winSize)
    deletedNames<-names(deletInd)
    strNames<-toString(deletedNames)
    lens<-lenSeqs[deletInd]
    strlens<-toString(lens)

    warning(paste("Sequences",strNames,"with length",strlens,"were deleted. Their lenghts were smaller than the window size"))

    if(length(label)==length(lenSeqs)){
      label<-label[-deletInd]
    }
    lenSeqs<-lenSeqs[-deletInd]
    seqs<-seqs[-deletInd]
  }

  dict<-list("A"=1,"C"=2,"D"=3,"E"=4,"F"=5,"G"=6,"H"=7,"I"=8,"K"=9,"L"=10,"M"=11,"N"=12,"P"=13,"Q"=14,"R"=15,"S"=16,"T"=17,"V"=18,"W"=19,"Y"=20)
  numSeqs<-length(seqs)


  if(overLap==FALSE){
    st<-0:(ceiling(lenSeqs[1]/winSize)-1)
    st<-(st*winSize)+1
    en<-st+winSize-1
    len=ceiling(lenSeqs[1]/winSize)
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


    featureMatrix<-matrix(0,ncol = (20*len),nrow = numSeqs)

    tname<-nameKmer(k=1,type = "aa")
    tname<-rep(tname,len)
    wname<-rep(0:(len-1),each =20)

    coln<-paste(tname,"w",wname,sep = "")

    colnames(featureMatrix)<-coln


    for(n in 1:numSeqs){

      seq<-seqs[n]
      subSeqs<-substring(seq,st,en)
      temp<-lapply(subSeqs, function(x){
        temp<-unlist(strsplit(x,""))
        tabtemp<-table(temp)
        tempvect<-vector(mode = "numeric",length = 20)
        names(tempvect)<-c("A","C","D" ,"E" ,"F" ,"G" ,"H" ,"I" ,"K" ,"L" ,"M" ,"N", "P" ,"Q" ,"R", "S" ,"T" ,"V" ,"W","Y")
        tempvect[names(tabtemp)]<-tabtemp
        return(tempvect)
      })
      temp2<-unlist(temp)
      featureMatrix[n,]<-temp2

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
      temp<-lapply(subSeqs, function(x){
        temp<-unlist(strsplit(x,""))
        tabtemp<-table(temp)
        tempvect<-vector(mode = "numeric",length = 20)
        names(tempvect)<-c("A","C","D" ,"E" ,"F" ,"G" ,"H" ,"I" ,"K" ,"L" ,"M" ,"N", "P" ,"Q" ,"R", "S" ,"T" ,"V" ,"W","Y")
        tempvect[names(tabtemp)]<-tabtemp
        return(tempvect)
      })
      temp2<-unlist(temp)
      names(temp2)<-NULL
      temp2<-c(nameSeq[n],temp2)
      temp3<-paste(temp2,collapse = "\t")
      write(temp3, outputFileDist, append=TRUE)
    }

  }
}


