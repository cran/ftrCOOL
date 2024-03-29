#' Di riboNucleotide To Binary RNA (DiNUC2Binary_RNA)
#'
#'
#' This function transforms a di-ribonucleotide to a binary number with four bits which is enough to represent all the possible types of di-ribonucleotides.
#' The type of the binary format is determined by the binaryType parameter.
#' For details about each format, please refer to the description of the binaryType parameter.
#'
#' @param binaryType It can take any of the following values: ('strBin','logicBin','numBin').
#' 'strBin' (String binary): each di-ribonucleotide is represented by a string containing 4 characters(0-1). For example, AA = "0000"   AC="0001"   ...  TT="1111"
#' 'logicBin' (logical value): Each di-ribonucleotide is represented by a vector containing 4 logical entries. For example, AA = c(F,F,F,F)   AC=c(F,F,F,T)  ...    TT=c(T,T,T,T)
#' 'numBin' (numeric bin): Each di-ribonucleotide is represented by a numeric (i.e., integer) vector containing 4 numeric entries. For example, AA = c(0,0,0,0)   AC = c(0,0,0,1)  ...  TT = c(1,1,1,1)
#'
#'
#' @note This function is provided for sequences with the same lengths.
#' Users can use 'txt' option in outFormat for sequences with different lengths.
#' Warning: If outFormat is set to 'mat' for sequences with different lengths, it returns an error.
#' Also, when output format is 'txt', label information is not shown in the text file.
#' It is noteworthy that 'txt' format is not usable for machine learning purposes if sequences have different sizes. Otherwise 'txt' format
#' is also usable for machine learning purposes.
#'
#'
#'
#' @param seqs is a FASTA file containing ribonucleotide sequences. The sequences start
#' with '>'. Also, seqs could be a string vector. Each element of the vector is a ribonucleotide sequence.
#'
#'
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @param outFormat (output format) can take two values: 'mat'(matrix) and 'txt'. The default value is 'mat'.
#'
#' @param outputFileDist shows the path and name of the 'txt' output file.
##'
#'
#' @return The output is different depending on the outFormat parameter ('mat' or 'txt').
#' If outFormat is set to 'mat', it returns a feature matrix for sequences with the same lengths.
#' The number of rows is equal to the number of sequences and if binaryType is 'strBin', the number of columns is (length of the sequences-1).
#' Otherwise, it is equal to (length of the sequences-1)*4.
#' If outFormat is 'txt', all binary values will be written to a 'txt' file. Each line in the file shows the binary format of a sequence.
#'
#'
#'
#' @export
#'
#' @examples
#'
#' fileLNC<-system.file("extdata/Carica_papaya101RNA.txt",package="ftrCOOL")
#' mat<-DiNUC2Binary_RNA(seqs = fileLNC, binaryType="numBin",outFormat="mat")
#'

DiNUC2Binary_RNA<-function(seqs,binaryType="numBin",outFormat="mat",outputFileDist="",label=c()){

  dict<-list("A"=1,"C"=2,"G"=3,"U"=4)
  nuc<-names(unlist(dict))

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

  numSeqs=length(seqs)

  lenSeqs<-sapply(seqs,nchar)




  featureVect<-vector()

  BinaryDiVect<-c("AA"="0000","AC"="0001","AG"="0010",
                  "AU"="0011","CA"="0100","CC"="0101",
                  "CG"="0110","CU"="0111","GA"="1000",
                  "GC"="1001","GG"="1010","GU"="1011",
                  "UA"="1100","UC"="1101","UG"="1110",
                  "UU"="1111")


  if(outFormat=="mat"){

    if(length(unique(lenSeqs))>1){
      stop("ERROR: All sequences should have the same length in 'mat' mode. For sequences with different lengths, please use 'txt' for outFormat parameter")
    }
    len<-lenSeqs[1]-1

    if(binaryType=="strBin"){
      featureMatrix<-matrix("",nrow = numSeqs,ncol = len)
      for(n in 1:numSeqs){
        seq<-seqs[n]
        chars<-unlist(strsplit(seq,""))

        temp1<-chars[1:len]
        temp2<-chars[2:(len+1)]
        dimer<-paste0(temp1,temp2)
        featureMatrix[n,]<-BinaryDiVect[dimer]
      }

      colnames(featureMatrix)<-paste("pos",1:len,sep = "")

    } else if(binaryType=="logicBin"){

      featureMatrix<-matrix(FALSE,nrow = numSeqs, ncol = (len*4))

      for(n in 1:numSeqs){
        seq<-seqs[n]
        chars<-unlist(strsplit(seq,split = ""))
        temp1<-chars[1:len]
        temp2<-chars[2:(len+1)]
        dimer<-paste0(temp1,temp2)
        bin<-BinaryDiVect[dimer]
        str<-paste(bin,collapse = '')
        strL<-unlist(strsplit(str,split = ''))
        strLogic<-as.logical(as.numeric(strL))
        featureMatrix[n,]<-strLogic
      }
      colnames(featureMatrix)<-paste("pos:",rep(1:len,each=4),"-",rep(nuc,len))

    } else if(binaryType=="numBin"){

      featureMatrix<-matrix(0,nrow = numSeqs, ncol = (len*4))

      for(n in 1:numSeqs){
        seq<-seqs[n]
        chars<-unlist(strsplit(seq,split = ""))
        temp1<-chars[1:len]
        temp2<-chars[2:(len+1)]
        dimer<-paste0(temp1,temp2)

        bin<-BinaryDiVect[dimer]
        str<-paste(bin,collapse = '')
        strL<-unlist(strsplit(str,split = ''))
        featureMatrix[n,]<-as.numeric(strL)

      }
      colnames(featureMatrix)<-paste("pos:",rep(1:len,each=4),"-",rep(nuc,len))

    } else{
      stop("ERROR! Choose one of 'strBin', 'logicBin', or 'numBin' for binaryFormat")
    }


    if(length(label)==numSeqs){
      featureMatrix<-as.data.frame(featureMatrix)
      featureMatrix<-cbind(featureMatrix,label)
    }
    row.names(featureMatrix)<-names(seqs)
    return(featureMatrix)

  } else if(outFormat=="txt"){
    nameSeq<-names(seqs)
    for(n in 1:numSeqs){
      seq<-seqs[n]
      chars<-unlist(strsplit(seq,""))
      temp1<-chars[1:(lenSeqs[n]-1)]
      temp2<-chars[2:lenSeqs[n]]
      dimer<-paste0(temp1,temp2)
      temp<-c(nameSeq[n],BinaryDiVect[dimer])
      temp<-paste(temp,collapse = "\t")
      write(temp,outputFileDist,append = TRUE)
    }

  }
  else{
    stop("ERROR: outFormat should be 'mat' or 'txt' ")
  }


}
