#' K Nucleotide Composition
#'
#' This function calculates the frequency of all k-mers in the sequence.
#'
#' @param seqs is a FASTA file containing nucleotide sequences. The sequences start
#' with '>'. Also, seqs could be a string vector. Each element of the vector is a nucleotide sequence.
#'
#'
#' @param rng This parameter can be a number or a vector. Each entry of the vector holds the value of k in the k-mer composition.
#' For each k in the rng vector, a new vector (whose size is 20^k) is created which contains the frequency of kmers.
#'
#' @param upto It is a logical parameter. The default value is FALSE. If rng is a number and upto is set to TRUE, rng is converted
#' to a vector with values from 1 to rng.
#'
#' @param reverse It is a logical parameter which assumes the reverse complement of the sequence.
#'
#' @param normalized is a logical parameter. When it is FALSE, the return value of the function does not change. Otherwise, the return value is normalized using the length of the sequence.
#'
#' @param ORF (Open Reading Frame) is a logical parameter. If it is set to true, ORF region of each sequence is considered instead of the original sequence (i.e., 3-frame).
#'
#' @param reverseORF is a logical parameter. It is enabled only if ORF is true.
#' If reverseORF is true, ORF region will be searched in the sequence and also in the reverse complement of the sequence (i.e., 6-frame).
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @return This function returns a feature matrix. The number of rows is equal to the number of sequences and
#' the number of columns depends on the rng vector. For each value k in the vector, (4)^k columns are created in the matrix.
#'
#' @export
#'
#' @examples
#'
#' fileLNC<-system.file("extdata/Athaliana_LNCRNA.fa",package="ftrCOOL")
#' mat<-kNUComposition(seqs=fileLNC,rng=c(1,3))
#'

kNUComposition<-function(seqs,rng=3,reverse=FALSE,upto=FALSE,normalized=TRUE,ORF=FALSE,reverseORF=TRUE,label=c()){


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

  if(ORF==TRUE){
    seqs=maxORF(seqs,reverse=reverseORF)
  }

  if(upto==TRUE && length(rng)==1){
    l<-length(rng)
    l<-rng[l]
    rng<-0:l
  }
  rng <- sort(rng)
  rng <- unique(rng)
  len<-length(rng)

  dict<-list("A"=1,"C"=2,"G"=3,"T"=4)
  mergedMatrix<-vector(mode = "numeric")
  numSeqs<-length(seqs)

  for(l in rng){

    featureMatrix<-matrix(0,ncol = (4^l),nrow = numSeqs)
    namesKmer<-nameKmer(l,type = "dna")
    colnames(featureMatrix)<-namesKmer

    if (l%%2==1)
    {
      size <- (4^l)/2
    } else{
      size <- ((4^l)/2)+(4^(l/2)/2)
    }

    for(n in 1:numSeqs){

      seq<-seqs[n]
      seqChars<-unlist(strsplit(seq,split = ""))
      lenSeq<-length(seqChars)
      kmers<-""

      #create all kmers occure in the seq
      for (i in 0:(l-1)){
        temp<-seqChars[(1+i):(lenSeq-(l-1-i))]
        kmers<-paste(kmers,temp,sep = "")
      }


      # table kmers of the seq
      tabKmers<-table(kmers)

      # a vector with name for each kmer
      tabNames<-names(tabKmers)

      #access to each kmer with convert its name to a number
      for(i in 1:length(tabKmers))
      {
        temp<-unlist(strsplit(tabNames[i],split = ""))
        num=0
        for(j in 1:l){
          pow<-4^(l-j)
          num<-num+(((as.numeric(dict[temp[j]]))-1)*pow)
        }
        num<-num+1
        featureMatrix[n,num]<-tabKmers[i]
      }

    }

    if(reverse==TRUE){

      for(i in (4^l):1){
        revCmpCl<-revComp(namesKmer[i])
        numRev=0
        for(j in 1:l){
          pow<-4^(l-j)
          numRev<-numRev+(((as.numeric(dict[revCmpCl[j]]))-1)*pow)
        }
        numRev<-numRev+1
        if(numRev<i){
          featureMatrix[,numRev]<-featureMatrix[,numRev]+featureMatrix[,i]
          featureMatrix<-featureMatrix[,-i]
        }
        if(ncol(featureMatrix)==size){
          break()
        }
      }

    }
    mergedMatrix<-cbind(mergedMatrix,featureMatrix)

  }

  if(normalized==TRUE){
    seqLen<-sapply(seqs, nchar)
    mergedMatrix<-mergedMatrix/seqLen
  }
  if(length(label)==numSeqs){
    mergedMatrix<-as.data.frame(mergedMatrix)
    mergedMatrix<-cbind(mergedMatrix,label)
  }
  row.names(mergedMatrix)<-names(seqs)


  return(mergedMatrix)
}


