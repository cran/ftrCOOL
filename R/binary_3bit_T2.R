#' Binary - 3bit - Type2 (binary_3bit_T2)
#'
#' This group of functions(binary_3bit_T1-T7) categorizes amino acids in 3 groups based on the type.
#' Then represent group of amino acids by a three dimentional vector.
#' The type of the binary format is determined by the binaryType parameter.
#' For details about each format, please refer to the description of the binaryType parameter.
#'
#' @param binaryType It can take any of the following values: ('strBin','logicBin','numBin').
#' 'strBin'(String binary): each amino acid is represented by a string containing 20 characters(0-1). For example, A = ALANIN = "1000000...0"
#' 'logicBin'(logical value): Each amino acid is represented by a vector containing 20 logical entries. For example, A = ALANIN = c(T,F,F,F,F,F,F,...F)
#' 'numBin' (numeric bin): Each amino acid is represented by a numeric (i.e., integer) vector containing 20 numerals. For example, A = ALANIN = c(1,0,0,0,0,0,0,...,0)
#'
#' @note This function is provided for sequences with the same lengths.
#' Users can use 'txt' option in outFormat for sequences with different lengths.
#' Warning: If outFormat is set to 'mat' for sequences with different lengths, it returns an error.
#' Also, when output format is 'txt', label information is not shown in the text file.
#' It is noteworthy that 'txt' format is not usable for machine learning purposes if sequences have different sizes. Otherwise 'txt' format
#' is also usable for machine learning purposes.
#'
#'
#' @param seqs is a FASTA file with amino acid sequences. Each sequence starts
#' with a '>' character. Also, seqs could be a string vector. Each element of the vector is a peptide/protein sequence.
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
#' @return The output is different depending on the outFormat parameter ('mat' or 'txt').
#' If outFormat is set to 'mat', it returns a feature matrix for sequences with the same lengths.
#' The number of rows is equal to the number of sequences and if binaryType is 'strBin', the number of columns is the length of the sequences.
#' Otherwise, it is equal to (length of the sequences)*3.
#' If outFormat is 'txt', all binary values will be written to a the output is written to a tab-delimited file. Each line in the file shows the binary format of a sequence.
#'
#' @export
#'
#'
#' @examples
#'
#' ptmSeqsADR<-system.file("extdata/",package="ftrCOOL")
#' ptmSeqsVect<-as.vector(read.csv(paste0(ptmSeqsADR,"/ptmVect101AA.csv"))[,2])
#' mat<-binary_3bit_T2(seqs = ptmSeqsVect, binaryType="numBin",outFormat="mat")
#'


binary_3bit_T2 <- function(seqs,binaryType="numBin",label=c(),outFormat="mat",outputFileDist="")
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
  numSeqs=length(seqs)


  Group=list("0-2.78"=c("G","A","S","T","P","D","C"),"2.95-4.0"=c("N","V","E","Q","I","L"),"4.03-8.08"=c("M","H","K","F","R","Y","W"))
  grps<-unlist(Group)

  dict<-vector()
  numGrp=3
  for (i in 1:numGrp)
  {
      vect<-rep(i,length(Group[[i]]))
      dict<-c(dict,vect)
  }
  names(dict)<-grps
  aa<-grps



  if(outFormat=="mat"){

    if(length(unique(lenSeqs))>1){
      stop("ERROR: All sequences should have the same length in 'mat' mode. For sequences with different lengths, please use 'txt' for outFormat parameter")
    }

    if(binaryType=="strBin"){

      binary<-rep(strrep(0,3),20)
      names(binary)=aa
      for(i in 1:length(dict))
      {
        pos<-dict[i]
        substr(binary[i],pos,pos)<-"1"
      }
      featureMatrix<-matrix("",nrow = numSeqs, ncol = lenSeqs[1])
      for(n in 1:numSeqs){
        seq<-seqs[n]
        charSeq<-unlist(strsplit(seq,split = ""))
        vect<-unlist(binary[charSeq])
        featureMatrix[n,]<-vect
      }

      colnames(featureMatrix)<-paste("pos",1:lenSeqs[1],sep = "")

    } else if(binaryType=="logicBin"){

      featureMatrix<-matrix(FALSE,nrow = numSeqs, ncol = (lenSeqs[1]*3))
      rng<-(0:(lenSeqs[1]-1))*3
      for(n in 1:numSeqs){
        seq<-seqs[n]
        charSeq<-unlist(strsplit(seq,split = ""))
        pos1<-as.numeric(dict[charSeq])
        pos1<-rng+pos1
        featureMatrix[n,pos1]<-TRUE
      }
      colnames(featureMatrix)<-paste0("pos",rep(1:lenSeqs[1],each=3),"-",rep(names(Group),lenSeqs[1]))

    } else if(binaryType=="numBin"){

      featureMatrix<-matrix(0,nrow = numSeqs, ncol = (lenSeqs[1]*3))
      rng<-(0:(lenSeqs[1]-1))*3
      for(n in 1:numSeqs){
        seq<-seqs[n]
        charSeq<-unlist(strsplit(seq,split = ""))
        pos1<-as.numeric(dict[charSeq])
        pos1<-rng+pos1
        featureMatrix[n,pos1]<-1
      }
      colnames(featureMatrix)<-paste0("pos",rep(1:lenSeqs[1],each=3),"-",rep(names(Group),lenSeqs[1]))

    } else{
      stop("ERROR! Choose one of 'strBin', 'logicBin', or 'numBin' for binaryFormat")
    }


    if(length(label)==numSeqs){
      featureMatrix<-as.data.frame(featureMatrix)
      featureMatrix<-cbind(featureMatrix,label)
    }
    row.names(featureMatrix)<-names(seqs)
    return(featureMatrix)
  }
  else if(outFormat=="txt"){

    vect<-vector()
    nameSeq<-names(seqs)
    binary<-rep(strrep(0,3),20)
    names(binary)=aa
    for(i in 1:length(dict))
    {
      pos<-dict[i]
      substr(binary[i],pos,pos)<-"1"
    }
    for(n in 1:numSeqs){
      seq<-seqs[n]
      charSeq<-unlist(strsplit(seq,split = ""))
      temp<-c(nameSeq[n],binary[charSeq])
      temp<-paste(temp,collapse = "\t")
      write(temp,outputFileDist,append = TRUE)
    }

  }
  else{
    stop("ERROR: outFormat should be 'mat' or 'txt' ")
  }

}

