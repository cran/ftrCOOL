#' Position-specific of four ribonucleotide (PS4_RNA)
#'
#' This function transforms each 4-ribonucleotide of the sequence to a binary format.
#' The type of the binary format is determined by the binaryType parameter.
#' For details about each format, please refer to the description of the binaryType parameter.
#'
#' @references Zhen Chen, Pei Zhao, Chen Li, Fuyi Li, Dongxu Xiang, Yong-Zi Chen, Tatsuya Akutsu, Roger J Daly, Geoffrey I Webb, Quanzhi Zhao, Lukasz Kurgan, Jiangning Song, iLearnPlus: a comprehensive and automated machine-learning platform for nucleic acid and protein sequence analysis, prediction and visualization, Nucleic Acids Research, (2021).
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
#' @param seqs is a FASTA file containing ribonucleotide sequences. The sequences start
#' with '>'. Also, seqs could be a string vector. Each element of the vector is a ribonucleotide sequence.
#'
#' @param binaryType It can take any of the following values: ('strBin','logicBin','numBin').
#' 'strBin'(String binary): each di-ribonucleotide is represented by a string containing 256 characters (255 times '0' and one '1'). For example, 'AAA' = "1000000000000000...0", ....
#' 'logicBin'(logical value): Each amino acid is represented by a vector containing 256 logical entries (255 times F and one T). For example, 'AA' = c(T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,...,F), ...
#' 'numBin' (numeric bin): Each amino acid is represented by a numeric (i.e., integer) vector containing 256 numerals (255 times '0' and one '1'). For example, 'AA' = c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,...,0)
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
#' Otherwise, it is equal to (length of the sequences-3)*256.
#' If outFormat is 'txt', all binary values will be written to a the output is written to a tab-delimited file. Each line in the file shows the binary format of a sequence.
#'
#' @export
#'
#'
#' @examples
#'
#' fileLNC<-system.file("extdata/Carica_papaya101RNA.txt",package="ftrCOOL")
#' mat<-PS4_RNA(seqs = fileLNC, binaryType="numBin",outFormat="mat")
#'
#'
#'
#'



PS4_RNA<- function(seqs,binaryType="numBin",label=c(),outFormat="mat",outputFileDist="")
{


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
  lenSeqs<-sapply(seqs,nchar)
  numSeqs=length(seqs)



  dict=1:256
  names(dict)=nameKmer(k=4,type="rna")

  nuc<-names(dict)
  if(outFormat=="mat"){

    if(length(unique(lenSeqs))>1){
      stop("ERROR: All sequences should have the same length in 'mat' mode. For sequences with different lengths, please use 'txt' for outFormat parameter")
    }
    lenSeq<-lenSeqs[1]

    if(binaryType=="strBin"){

      binary<-rep(strrep(0,256),256)
      names(binary)=names(dict)
      for(i in 1:length(dict))
      {
        pos<-dict[i]
        substr(binary[i],pos,pos)<-"1"
      }
      featureMatrix<-matrix("",nrow = numSeqs, ncol = (lenSeq-3))
      for(n in 1:numSeqs){
        seq<-seqs[n]
        seqChars<-unlist(strsplit(seq,split = ""))
        temp<-seqChars[1:(lenSeq-3)]
        temp2<-seqChars[2:(lenSeq-2)]
        temp3<-seqChars[3:(lenSeq-1)]
        temp4<-seqChars[4:lenSeq]
        Fmer<-paste(temp,temp2,temp3,temp4,sep = "")
        vect<-unlist(binary[Fmer])
        featureMatrix[n,]<-vect
      }

      colnames(featureMatrix)<-paste("pos",1:(lenSeq-3),sep = "")

    } else if(binaryType=="logicBin"){

      featureMatrix<-matrix(FALSE,nrow = numSeqs, ncol = ((lenSeq-3)*256))
      rng<-(0:(lenSeq-4))*256
      for(n in 1:numSeqs){
        seq<-seqs[n]

        seqChars<-unlist(strsplit(seq,split = ""))
        temp<-seqChars[1:(lenSeq-3)]
        temp2<-seqChars[2:(lenSeq-2)]
        temp3<-seqChars[3:(lenSeq-1)]
        temp4<-seqChars[4:lenSeq]
        Fmer<-paste(temp,temp2,temp3,temp4,sep = "")


        pos1<-as.numeric(dict[Fmer])
        pos1<-rng+pos1
        featureMatrix[n,pos1]<-TRUE
      }
      colnames(featureMatrix)<-paste("pos:",rep(1:(lenSeq-3),each=256),"-",rep(nuc,(lenSeq-3)))

    } else if(binaryType=="numBin"){

      featureMatrix<-matrix(0,nrow = numSeqs, ncol = ((lenSeq-3)*256))
      rng<-(0:(lenSeq-4))*256
      for(n in 1:numSeqs){
        seq<-seqs[n]

        seqChars<-unlist(strsplit(seq,split = ""))
        temp<-seqChars[1:(lenSeq-3)]
        temp2<-seqChars[2:(lenSeq-2)]
        temp3<-seqChars[3:(lenSeq-1)]
        temp4<-seqChars[4:lenSeq]
        Fmer<-paste(temp,temp2,temp3,temp4,sep = "")

        pos1<-as.numeric(dict[Fmer])
        pos1<-rng+pos1
        featureMatrix[n,pos1]<-1
      }
      colnames(featureMatrix)<-paste("pos:",rep(1:(lenSeq-3),each=256),"-",rep(nuc,(lenSeq-3)))

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
    binary<-rep(strrep(0,256),256)
    names(binary)=nuc
    for(i in 1:length(dict))
    {
      pos<-dict[i]
      substr(binary[i],pos,pos)<-"1"
    }
    for(n in 1:numSeqs){
      seq<-seqs[n]

      seqChars<-unlist(strsplit(seq,split = ""))
      temp<-seqChars[1:(lenSeq-3)]
      temp2<-seqChars[2:(lenSeq-2)]
      temp3<-seqChars[3:(lenSeq-1)]
      temp4<-seqChars[4:lenSeq]
      Fmer<-paste(temp,temp2,temp3,temp4,sep = "")

      temp<-c(nameSeq[n],binary[Fmer])
      temp<-paste(temp,collapse = "\t")
      write(temp,outputFileDist,append = TRUE)
    }

  }
  else{
    stop("ERROR: outFormat should be 'mat' or 'txt' ")
  }

}

