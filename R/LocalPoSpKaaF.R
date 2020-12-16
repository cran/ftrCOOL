#' Local Position Specific k Amino Acids Frequency (LocalPoSpKAAF)
#'
#' For each sequence, this function creates a feature vector denoted as
#' (f1,f2, f3, …, fN), where fi = freq(i'th k-mer of the sequence) / i.
#' It should be applied to sequences with the same length.
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
#' @param k is a numeric value which holds the value of k in the k-mers.
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @param outFormat (output format) can take two values: 'mat'(matrix) and 'txt'. The default value is 'mat'.
#'
#' @param outputFileDist shows the path and name of the 'txt' output file.
#'
#' @return The output depends on the outFormat parameter which can be either 'mat' or 'txt'. If outFormat is 'mat', the function returns a feature
#' matrix for sequences with the same length such that the number of columns is (sequence length-k+1)
#' and the number of rows is equal to the number of sequences.
#' If the outFormat is 'txt', the output is written to a tab-delimited file.
#'
#' @export
#'
#' @examples
#'
#' dir = tempdir()
#' ptmSeqsADR<-system.file("extdata/",package="ftrCOOL")
#' ptmSeqsVect<-as.vector(read.csv(paste0(ptmSeqsADR,"/ptmVect101AA.csv"))[,2])
#' mat<-LocalPoSpKAAF(seqs = ptmSeqsVect, k=2,outFormat="mat")
#'
#' ad<-paste0(dir,"/LocalPoSpKaaF.txt")
#' filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
#' LocalPoSpKAAF(seqs = filePrs, k=1,outFormat="txt"
#' ,outputFileDist=ad)
#'

LocalPoSpKAAF<-function(seqs,k=2,label=c(),outFormat="mat",outputFileDist=""){


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
  numSeqs<-length(seqs)
  lenSeqs<-sapply(seqs, nchar)


  if(outFormat=="mat"){
  if(length(unique(lenSeqs))>1){
    stop("ERROR: All sequences should have the same length in 'mat' mode. For sequences with different lengths, please use 'txt' for outFormat parameter")
  }

  featureMatrix<-matrix(0,nrow = numSeqs,ncol = (lenSeqs[1]-k+1))
  colnames(featureMatrix)<-paste0("f",rep(1:ncol(featureMatrix)))


  for(n in 1:numSeqs){

    seq<-seqs[n]
    seqChars<-unlist(strsplit(seq,split = ""))


    numVect<-vector()
    #creat kmers occure in the sequence
    kmers<-""
    for (i in 0:(k-1)){
      temp<-seqChars[(1+i):(lenSeqs[n]-(k-1-i))]
      kmers<-paste(kmers,temp,sep = "")
    }

    for(i in 1:length(kmers)){

      subseq<-kmers[i]
      a=numVect[subseq]
      if(is.na(a)){
        numVect<-c(numVect,1)
        names(numVect)[length(numVect)]<-subseq
        featureMatrix[n,i]<-numVect[subseq]/(i+k-1)

      } else {
        numVect[subseq]<-numVect[subseq]+1
        featureMatrix[n,i]<-numVect[subseq]/(i+k-1)

      }

    }

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
      ve<-vector(mode = "numeric",length = lenSeqs[n])
      seq<-seqs[n]
      seqChars<-unlist(strsplit(seq,split = ""))


      numVect<-vector()
      #creat kmers occure in the sequence
      kmers<-""
      for (i in 0:(k-1)){
        temp<-seqChars[(1+i):(lenSeqs[n]-(k-1-i))]
        kmers<-paste(kmers,temp,sep = "")
      }
      lenkmer<-length(kmers)
      vect<-vector(mode = "numeric",length = (lenSeqs[n]-k+1))
      for(i in 1:length(kmers)){
        #subseq<-substr(seq,i,(i+k-1))
        subseq<-kmers[i]
        a=numVect[subseq]
        if(is.na(a)){
          numVect<-c(numVect,1)
          names(numVect)[length(numVect)]<-subseq
          vect[i]<-numVect[subseq]/(i+k-1)

        } else {
          numVect[subseq]<-numVect[subseq]+1
          vect[i]<-numVect[subseq]/(i+k-1)

        }

      }
      vect<-c(nameSeq[n],vect)
      vect2<-paste(vect,collapse = "\t")
      write(vect2,outputFileDist,append = TRUE)

    }


  }
  else {
    stop("ERROR: outFormat should be 'mat' or 'txt' ")
  }

}
