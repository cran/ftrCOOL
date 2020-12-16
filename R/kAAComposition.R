#' k Amino Acid Composition (kAAComposition)
#'
#' This function calculates the frequency of all k-mers in the sequence(s).
#'
#'
#' @param seqs is a FASTA file with amino acid sequences. Each sequence starts
#' with a '>' character. Also, seqs could be a string vector. Each element of the vector is a peptide/protein sequence.
#'
#'
#' @param rng This parameter can be a number or a vector. Each entry of the vector holds the value of k in the k-mer composition.
#'
#' @param upto It is a logical parameter. The default value is FALSE. If rng is a number and upto is set to TRUE, rng is converted
#' to a vector with values from 1 to rng.
#'
#' @param normalized is a logical parameter. When it is FALSE, the return value of the function does not change. Otherwise, the return value is normalized using the length of the sequence.
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @return This function returns a feature matrix. The number of rows is equal to the number of sequences and
#' the number of columns depends on rng vector. For each value k in the vector, (20)^k columns are created in the matrix.
#'
#'
#' @export
#'
#' @examples
#'
#' filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
#'
#' mat1<-kAAComposition(seqs=filePrs,rng=3,upto=TRUE)
#' mat2<-kAAComposition(seqs=filePrs,rng=c(1,3),upto=TRUE)


kAAComposition<-function(seqs,rng=3,upto=FALSE,normalized=TRUE,label=c()){

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

  if(upto==TRUE && length(rng)==1){
    l<-length(rng)
    l<-rng[l]
    rng<-0:l
  }
  rng <- sort(rng)
  rng <- unique(rng)
  len<-length(rng)

  dict<-list("A"=1,"C"=2,"D"=3,"E"=4,"F"=5,"G"=6,"H"=7,"I"=8,"K"=9,"L"=10,"M"=11,"N"=12,"P"=13,"Q"=14,"R"=15,"S"=16,"T"=17,"V"=18,"W"=19,"Y"=20)
  mergedMatrix<-vector(mode = "numeric")
  numSeqs<-length(seqs)

  for(l in rng){

    featureMatrix<-matrix(0,ncol = (20^l),nrow = numSeqs)
    namesKmer<-nameKmer(l,type = "aa")
    colnames(featureMatrix)<-namesKmer

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
          pow<-20^(l-j)
          num<-num+(((as.numeric(dict[temp[j]]))-1)*pow)
        }
        num<-num+1
        featureMatrix[n,num]<-tabKmers[i]
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


