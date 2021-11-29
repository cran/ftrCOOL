#' Binary - 6bit (binary_6bit)
#'
#' This function categorizes amino acids in 6 groups.
#' Then represent group of amino acids by a 6 dimentional vector i.e.{e1, e2, e3, e4, e5, e6}.
#' e1={H, R, K}, e2={D, E, N, D}, e3={C}, e4={S, T, P, A, G}, e5={M, I, L, V}, e6={F, Y, W}.
#' e1 is ecoded by 100000 e2 is encoded by 010000 and ... and e6 is encoded by 000001.
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
#' Otherwise, it is equal to (length of the sequences)*6.
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
#' mat<-binary_6bit(seqs = ptmSeqsVect, binaryType="numBin",outFormat="mat")
#'


binary_6bit <- function(seqs,binaryType="numBin",label=c(),outFormat="mat",outputFileDist="")
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
  #dict<-list("A"=1,"C"=2,"D"=3,"E"=4,"F"=5,"G"=6,"H"=7,"I"=8,"K"=9,"L"=10,"M"=11,"N"=12,"P"=13,"Q"=14,"R"=15,"S"=16,"T"=17,"V"=18,"W"=19,"Y"=20)
  #aa<-names(unlist(dict))


  Group=list("Grp1"=c("H","R","K"),"Grp2"=c("Q","E","N","D"),"Grp3"=c("C"),"Grp4"=c("S","T","P","A","G"),"Grp5"=c("M","I","L","V"),"Grp6"=c("F","Y","W"))
  grps<-unlist(Group)

  dict<-vector()
  numGrp=6
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

      binary<-rep(strrep(0,6),20)
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

      featureMatrix<-matrix(FALSE,nrow = numSeqs, ncol = (lenSeqs[1]*6))
      rng<-(0:(lenSeqs[1]-1))*6
      for(n in 1:numSeqs){
        seq<-seqs[n]
        charSeq<-unlist(strsplit(seq,split = ""))
        pos1<-as.numeric(dict[charSeq])
        pos1<-rng+pos1
        featureMatrix[n,pos1]<-TRUE
      }
      #colnames(featureMatrix)<-paste("pos:",rep(1:lenSeqs[1],each=20),"-",rep(aa,lenSeqs[1]))

    } else if(binaryType=="numBin"){

      featureMatrix<-matrix(0,nrow = numSeqs, ncol = (lenSeqs[1]*6))
      rng<-(0:(lenSeqs[1]-1))*6
      for(n in 1:numSeqs){
        seq<-seqs[n]
        charSeq<-unlist(strsplit(seq,split = ""))
        pos1<-as.numeric(dict[charSeq])
        pos1<-rng+pos1
        featureMatrix[n,pos1]<-1
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
    binary<-rep(strrep(0,6),20)
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

