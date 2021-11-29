#' Nucleotide To Binary (NUC2Binary_DNA)
#'
#'
#' This function transforms a nucleotide to a binary format.
#' The type of the binary format is determined by the binaryType parameter.
#' For details about each format, please refer to the description of the binaryType parameter.
#'
#'
#' @note This function is provided for sequences with the same lengths.
#' Users can use 'txt' option in outFormat parameter for sequences with different lengths.
#' Warning: If outFormat is set to 'mat' for sequences with different lengths, it returns an error.
#' Also, when output format is 'txt', label information is not shown in the text file.
#' It is noteworthy that 'txt' format is not usable for machine learning purposes.
#'
#'
#' @param seqs is a FASTA file containing nucleotide sequences. The sequences start
#' with '>'. Also, seqs could be a string vector. Each element of the vector is a nucleotide sequence.
#'
#' @param binaryType It can take any of the following values: ('strBin','logicBin','numBin').
#' 'strBin'(String binary): each nucleotide is represented by a string containing 4 characters(0-1). A = "0001" , C = "0010" , G = "0100" , T = "1000"
#' 'logicBin'(logical value): Each nucleotide is represented by a vector containing 4 logical entries. A = c(F,F,F,T) , ... , T = c(T,F,F,F)
#' 'numBin' (numeric bin): Each nucleotide is represented by a numeric (i.e., integer) vector containing 4 numerals. A = c(0,0,0,1) , ... , T = c(1,0,0,0)
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
#' Otherwise, it is equal to (length of the sequences)*4.
#' If outFormat is 'txt', all binary values will be written to a 'txt' file. Each line in the file shows the binary format of a sequence.
#'
#'
#' @export
#'
#' @examples
#'
#' dir = tempdir()
#' LNCSeqsADR<-system.file("extdata/",package="ftrCOOL")
#' LNC50Nuc<-as.vector(read.csv(paste0(LNCSeqsADR,"/LNC50Nuc.csv"))[,2])
#' mat<-NUC2Binary_DNA(seqs = LNC50Nuc,outFormat="mat")
#'
#' ad<-paste0(dir,"/NUC2Binary.txt")
#' fileLNC<-system.file("extdata/Athaliana_LNCRNA.fa",package="ftrCOOL")
#' NUC2Binary_DNA(seqs = fileLNC,binaryType="numBin",outFormat="txt",outputFileDist=ad)
#'
#' unlink("dir", recursive = TRUE)


NUC2Binary_DNA <- function(seqs,binaryType="numBin",label=c(),outFormat="mat",outputFileDist="")
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
  numSeqs=length(seqs)
  dict<-list("A"=1,"C"=2,"G"=3,"T"=4)
  nuc=names(unlist(dict))


  if(outFormat=="mat"){

  if(length(unique(lenSeqs))>1){
    stop("ERROR: All sequences should have the same length in 'mat' mode. For sequences with different lengths, please use 'txt' for outFormat parameter")
  }

    if(binaryType=="strBin"){
      featureMatrix<-matrix(0,nrow = numSeqs, ncol = lenSeqs[1])
      colnames(featureMatrix)<-paste("pos",1:lenSeqs[1],sep = "")
      for(n in 1:numSeqs){
        chars<-unlist(strsplit(seqs[n],""))
        chars[chars=="A"]="0001"
        chars[chars=="C"]="0010"
        chars[chars=="G"]="0100"
        chars[chars=="T"]="1000"
        featureMatrix[n,]=chars
      }

    } else if(binaryType=="logicBin"){

      featureMatrix<-matrix(FALSE,nrow = numSeqs, ncol = (lenSeqs[1]*4))
      rng<-(0:(lenSeqs[1]-1))*4
      for(n in 1:numSeqs){
        seq<-seqs[n]
        charSeq<-unlist(strsplit(seq,split = ""))
        pos1<-as.numeric(dict[charSeq])
        pos1<-rng+pos1
        featureMatrix[n,pos1]<-TRUE
      }
      colnames(featureMatrix)<-paste("pos:",rep(1:lenSeqs[1],each=4),"-",rep(nuc,lenSeqs[1]))

    } else if(binaryType=="numBin"){

      featureMatrix<-matrix(0,nrow = numSeqs, ncol = (lenSeqs[1]*4))
      rng<-(0:(lenSeqs[1]-1))*4
      for(n in 1:numSeqs){
        seq<-seqs[n]
        charSeq<-unlist(strsplit(seq,split = ""))
        pos1<-as.numeric(dict[charSeq])
        pos1<-rng+pos1
        featureMatrix[n,pos1]<-1
      }
      colnames(featureMatrix)<-paste("pos:",rep(1:lenSeqs[1],each=4),"-",rep(nuc,lenSeqs[1]))

    } else{
      stop("ERROR: Choose one of 'strBin', 'logicBin', or 'numBin' for binaryFormat")
    }


  if(length(label)==numSeqs){
    featureMatrix<-as.data.frame(featureMatrix)
    featureMatrix<-cbind(featureMatrix,label)
  }
  row.names(featureMatrix)<-names(seqs)



  return(featureMatrix)
  }   else if(outFormat=="txt"){

    nameSeq<-names(seqs)
    for(n in 1:numSeqs){

      chars<-unlist(strsplit(seqs[n],""))
      chars[chars=="A"]="0001"
      chars[chars=="C"]="0010"
      chars[chars=="G"]="0100"
      chars[chars=="T"]="1000"

      temp<-c(nameSeq[n],chars)
      temp<-paste(temp,collapse = "\t")
      write(temp,outputFileDist,append = TRUE)
    }


  }
  else {
    stop("ERROR: outFormat should be 'mat' or 'txt' ")
  }


}

