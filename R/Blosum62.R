#' Blosum62 (BLOSUM62)
#'
#' This function creates a 20-dimentional numeric vector for each amino acid of a sequence. Each entry of the vector contains
#' the similarity score of the amino acid with other amino acids including itself. The score is extracted from the Blosum62 matrix.
#'
#'
#' @note This function is provided for sequences with the same lengths.
#' Users can use 'txt' option in outFormat for sequences with different lengths.
#' Warning: If outFormat is set to 'mat' for sequences with different lengths, it returns an error.
#' Also, when output format is 'txt', label information is not shown in the text file.
#' It is noteworthy that 'txt' format is not usable for machine learning purposes if sequences have different sizes. Otherwise 'txt' format
#' is also usable for machine learning purposes.
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
#' @return The output depends on the outFormat parameter which can be either 'mat' or 'txt'. If outFormat is 'mat', the function returns a feature
#' matrix for sequences with the same length such that the number of columns is (sequence length)*20
#' and the number of rows is equal to the number of sequences.
#' If the outFormat is 'txt', the output is written to a tab-delimited file.
#'
#'
#' @export
#'
#' @examples
#'
#' dir = tempdir()
#' ptmSeqsADR<-system.file("extdata/",package="ftrCOOL")
#' filePr<-system.file("extdata/protein.fasta",package="ftrCOOL")
#' filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
#'
#' ad<-paste0(dir,"/blosum62.txt")
#' vect<-BLOSUM62(seqs = filePr,outFormat="mat")
#' BLOSUM62(seqs = filePrs,outFormat="txt",outputFileDist=ad)
#'
#' unlink("dir", recursive = TRUE)
#'
#'
BLOSUM62<-function(seqs,label=c(),outFormat="mat",outputFileDist=""){


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



  featureMatrix<-matrix(0, nrow = numSeqs, ncol = (lenSeqs[1]*23))


  adBlosum<-paste0(path.pack,"/BLOSUM62.csv")
  Blosum<-read.csv(adBlosum)
  row.names(Blosum)<-Blosum[,1]
  Blosum<-Blosum[,-1]
  Blosum<-as.matrix(Blosum)
  Blosum<-type.convert(Blosum)
  nameP1<-rep((1:(lenSeqs[1])),each=23)
  nameP1<-paste0("aa",nameP1)
  aaVect<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","B","Z","X")

  if(outFormat=="mat"){

    if(length(unique(lenSeqs))>1){
      stop("ERROR: All sequences should have the same length in 'mat' mode. For sequences with different lengths, please use 'txt' for outFormat parameter")
    }
    nameP2<-paste("d",rep(aaVect,lenSeqs[1]),sep = "")
    colnames(featureMatrix)<-paste(nameP1,nameP2)
    for(n in 1:numSeqs){
      seq=seqs[n]
      chars=unlist(strsplit(seq,""))
      featureVector<-Blosum[,chars]
      featureVector<-as.vector(featureVector)
      featureMatrix[n,]<-featureVector
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
      featureVector<-Blosum[,chars]
      vect<-as.vector(featureVector)
      temp<-c(nameSeq[n],vect)
      temp<-paste(temp,collapse = "\t")
      write(temp,outputFileDist,append = TRUE)
    }
  }

}



