#' Enhanced Nucleotide Composition (ENUComposition_DNA)
#'
#' This function slides a window over the input sequence(s).
#' Also, it computes the composition of nucleotides that appears within the limits of the window.
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
#' and the number of rows is equal to the number of sequences.
#' If the outFormat is 'txt', the output is written to a tab-delimited file.
#'
#' @export
#'
#' @examples
#'
#' dir = tempdir()
#' LNCSeqsADR<-system.file("extdata/",package="ftrCOOL")
#' LNC50Nuc<-as.vector(read.csv(paste0(LNCSeqsADR,"/LNC50Nuc.csv"))[,2])
#' mat<-ENUComposition_DNA(seqs = LNC50Nuc, winSize=20,outFormat="mat")
#'
#' ad<-paste0(dir,"/ENUCcompos.txt")
#' fileLNC<-system.file("extdata/Athaliana_LNCRNA.fa",package="ftrCOOL")
#' ENUComposition_DNA(seqs = fileLNC,outFormat="txt",winSize=20
#' ,outputFileDist=ad,overLap=FALSE)


ENUComposition_DNA <- function(seqs,winSize=50,overLap=TRUE,label=c(),outFormat='mat',outputFileDist="")
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

    warning(paste("Sequences",strNames,"with length",strlens,"were deleted. Their lenghts were smaller than the window size"))

    if(length(label)==length(lenSeqs)){
      label<-label[-deletInd]
    }
    lenSeqs<-lenSeqs[-deletInd]
    seqs<-seqs[-deletInd]
  }

  numSeqs<-length(seqs)

  if(overLap==FALSE){
    st<-0:(ceiling(lenSeqs[1]/winSize)-1)
    st<-(st*winSize)+1
    en<-st+winSize-1
    len=ceiling(lenSeqs[1]/winSize)
  }
  else{
    st<-1:(lenSeqs[1]-winSize+1)
    en<-winSize:(lenSeqs[1])
    len<-lenSeqs[1]-winSize+1
  }

  if(outFormat=='mat'){
    if(length(unique(lenSeqs))>1){
      stop("ERROR: All sequences should have the same length in 'mat' mode. For sequences with different lengths, please use 'txt' for outFormat parameter")
    }


    featureMatrix<-matrix(0,ncol = (4*len),nrow = numSeqs)

    tname<-nameKmer(k=1,type = "dna")
    tname<-rep(tname,len)
    wname<-rep(0:(len-1),each =4)

    coln<-paste(tname,"w",wname,sep = "")

    colnames(featureMatrix)<-coln

    tempvect<-vector(mode = "numeric",length = 4)
    names(tempvect)<-c("A","C","G" ,"T")
    subSeqs<-lapply(seqs, function(x){
      subSeqs<-substring(x,st,en)
      return(subSeqs)
    })
    tempvect<-vector(mode = "numeric",length = 4)
    names(tempvect)<-c("A","C","G","T")
    listvect<-list()

    if(overLap){

      for(j in 1:numSeqs) {
        totalvect<-c()
        chars<-unlist(strsplit(subSeqs[[j]][1],""))
        for(i in 1:winSize){
          tempvect[chars[i]]<-tempvect[chars[i]]+1
        }
        totalvect<-c(totalvect,tempvect)
        for(i in 2:len){
          firstChar=substr(subSeqs[[j]][(i-1)],1,1)
          lastChar=substr(subSeqs[[j]][i],winSize,winSize)
          tempvect[firstChar]=tempvect[firstChar]-1
          tempvect[lastChar]=tempvect[lastChar]+1
          totalvect<-c(totalvect,tempvect)
        }
        featureMatrix[j,]<-totalvect
      }


    } else{

      totalvect<-c()
      for(j in 1:numSeqs) {
        tabsubs<-lapply(subSeqs[[j]], function(x){
          tempvect[]<-c(0,0,0,0)
          chars<-unlist(strsplit(x,""))
          for(i in 1:length(chars)){
            tempvect[chars[i]]=tempvect[chars[i]]+1
          }
          totalvect<-c(totalvect,tempvect)
          return(totalvect)
        })
        #listvect[[j]]<-unlist(tabsubs)
        featureMatrix[j,]<-unlist(tabsubs)
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
    totalvect<-c()
    nameSeq<-names(seqs)

    tempvect<-vector(mode = "numeric",length = 4)
    names(tempvect)<-c("A","C","G" ,"T")
    subSeqs<-lapply(seqs, function(x){
      subSeqs<-substring(x,st,en)
      return(subSeqs)
    })

    for(j in 1:numSeqs) {
      tabsubs<-lapply(subSeqs[[j]], function(x){
        tempvect[]<-c(0,0,0,0)
        chars<-unlist(strsplit(x,""))
        for(i in 1:length(chars)){
          tempvect[chars[i]]=tempvect[chars[i]]+1
        }
        totalvect<-c(totalvect,tempvect)
        return(totalvect)
      })
      temp2<-unlist(tabsubs)
      names(temp2)<-NULL
      temp2<-c(nameSeq[j],temp2)
      temp3<-paste(temp2,collapse = "\t")
      write(temp3, outputFileDist, append=TRUE)
    }

  }
}


