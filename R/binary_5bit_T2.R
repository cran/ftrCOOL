#' Binary - 5bit - Type2 (binary_5bit_T2)
#'
#' The idea behind this function is: We have 20 amino acids and we can show them with at least 5 bits.
#' A is encoded by (00011), C (00101), D (00110), E (00111), F(01001), G (01010), H (01011), I (01100), K (01101), L (01110), M (10001), N (10010), P (10011), Q (10100), R (10101), S (10110), T (11000), V (11001), W (11010), Y (11100).
#' This function transforms an amino acid to a binary format.
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
#' Otherwise, it is equal to (length of the sequences)*5.
#' If outFormat is 'txt', all binary values will be written to a the output is written to a tab-delimited file. Each line in the file shows the binary format of a sequence.
#'
#' @export
#'
#'
#'
#' @examples
#'
#' ptmSeqsADR<-system.file("extdata/",package="ftrCOOL")
#' ptmSeqsVect<-as.vector(read.csv(paste0(ptmSeqsADR,"/ptmVect101AA.csv"))[,2])
#' mat<-binary_5bit_T2(seqs = ptmSeqsVect, binaryType="numBin",outFormat="mat")
#'


binary_5bit_T2 <- function(seqs,binaryType="numBin",label=c(),outFormat="mat",outputFileDist="")
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

  dict<-list("A"=c(0,0,0,1,1),"C"=c(0,0,1,0,1),"D"=c(0,0,1,1,0),"E"=c(0,0,1,1,1),"F"=c(0,1,0,0,1),"G"=c(0,1,0,1,0),"H"=c(0,1,0,1,1),"I"=c(0,1,1,0,0),"K"=c(0,1,1,0,1),"L"=c(0,1,1,1,0),"M"=c(1,0,0,0,1),"N"=c(1,0,0,1,0),"P"=c(1,0,0,1,1),"Q"=c(1,0,1,0,0),"R"=c(1,0,1,0,1),"S"=c(1,0,1,1,0),"T"=c(1,1,0,0,0),"V"=c(1,1,0,0,1),"W"=c(1,1,0,1,0),"Y"=c(1,1,1,0,0))

  aa<-names(dict)

  if(outFormat=="mat"){

    if(length(unique(lenSeqs))>1){
      stop("ERROR: All sequences should have the same length in 'mat' mode. For sequences with different lengths, please use 'txt' for outFormat parameter")
    }

    if(binaryType=="strBin"){
      newdict<-unlist(lapply(dict,function(x){paste0(x,collapse = "")}))

      featureMatrix<-matrix("",nrow = numSeqs, ncol = lenSeqs[1])
      for(n in 1:numSeqs){
        seq<-seqs[n]
        charSeq<-unlist(strsplit(seq,split = ""))
        vect<-newdict[charSeq]
        featureMatrix[n,]<-vect
      }

      colnames(featureMatrix)<-paste("pos",1:lenSeqs[1],sep = "")

    } else if(binaryType=="logicBin"){

      featureMatrix<-matrix(FALSE,nrow = numSeqs, ncol = (lenSeqs[1]*5))

      for(n in 1:numSeqs){
        seq<-seqs[n]
        charSeq<-unlist(strsplit(seq,split = ""))
        vect<-unlist(dict[charSeq])
        vect<-as.logical(vect)
        featureMatrix[n,]<-vect

      }
      #colnames(featureMatrix)<-paste("pos:",rep(1:lenSeqs[1],each=20),"-",rep(aa,lenSeqs[1]))

    } else if(binaryType=="numBin"){

      featureMatrix<-matrix(0,nrow = numSeqs, ncol = (lenSeqs[1]*5))
      for(n in 1:numSeqs){
        seq<-seqs[n]
        charSeq<-unlist(strsplit(seq,split = ""))
        vect<-unlist(dict[charSeq])
        featureMatrix[n,]<-vect
      }
      #colnames(featureMatrix)<-paste("pos:",rep(1:lenSeqs[1],each=20),"-",rep(aa,lenSeqs[1]))

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
    newdict<-unlist(lapply(dict,function(x){paste0(x,collapse = "")}))
    for(n in 1:numSeqs){
      seq<-seqs[n]
      charSeq<-unlist(strsplit(seq,split = ""))
      temp<-c(nameSeq[n],newdict[charSeq])
      temp<-paste(temp,collapse = "\t")
      write(temp,outputFileDist,append = TRUE)
    }

  }
  else{
    stop("ERROR: outFormat should be 'mat' or 'txt' ")
  }

}

