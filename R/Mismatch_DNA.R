#' Mismatch_DNA (Mismatch_DNA)
#'
#' This function also calculates the frequencies of all k-mers in the sequence but alows maximum m mismatch.
#' m<k.
#'
#' @references Liu, B., Gao, X. and Zhang, H. BioSeq-Analysis2.0: an updated platform for analyzing DNA, RNA and protein sequences at sequence level and residue level based on machine learning approaches. Nucleic Acids Res (2019).
#'
#' @param seqs is a FASTA file containing nucleotide sequences. The sequences start
#' with '>'. Also, seqs could be a string vector. Each element of the vector is a nucleotide sequence.
#'
#'
#' @param k This parameter can be a number which shows kmer.
#'
#' @param m This parametr shows muximum number of mismatches.
#'
#'
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
#' mat<-Mismatch_DNA(seqs=fileLNC)
#'
#'

Mismatch_DNA<-function(seqs,k=3,m=2,label=c()){


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

  dict<-list("A"=1,"C"=2,"G"=3,"T"=4)
  mergedMatrix<-vector(mode = "numeric")
  numSeqs<-length(seqs)


    featureMatrix<-matrix(0,ncol = (4^k),nrow = numSeqs)
    namesKmer<-nameKmer(k,type = "dna")
    colnames(featureMatrix)<-namesKmer


    for(n in 1:numSeqs){

      seq<-seqs[n]
      seqChars<-unlist(strsplit(seq,split = ""))
      lenSeq<-length(seqChars)
      kmers<-""

      #create all kmers occure in the seq
      for (i in 0:(k-1)){
        temp<-seqChars[(1+i):(lenSeq-(k-1-i))]
        kmers<-paste(kmers,temp,sep = "")
      }


      # table kmers of the seq
      tabKmers<-table(kmers)

      # a vector with name for each kmer
      tabNames<-names(tabKmers)

      featureMatrix[n,tabNames]<-tabKmers


    }

    listChars<-lapply(namesKmer,function(x) {unlist(strsplit(x,""))})
    listMis<-list()
    len<-length(listChars)
    for(i in 1:len){
      newList<-listChars[-i]
      vect<-vector()
      for(j in 1:(len-1)){
        bltemp<-newList[[j]]==listChars[[i]]
        tabbltemp<-table(bltemp)
        if(tabbltemp["FALSE"]<=m){
          vect<-c(vect,paste(newList[[j]],collapse = ""))
        }
      }

      listMis[[i]]<-vect
    }
    mismatchMat<-matrix(0,ncol = 4^k,nrow = numSeqs)
    colnames(mismatchMat)<-colnames(featureMatrix)
    rownames(mismatchMat)<-names(seqs)
    if(nrow(featureMatrix)>1){
      for(i in 1:len){
        sumvect<-apply(featureMatrix[,listMis[[i]]], 1, sum)
        mismatchMat[,i]<-featureMatrix[,i]+sumvect
      }
    }
    else{
      for(i in 1:len){
        sumvect<-sum(featureMatrix[,listMis[[i]]])
        mismatchMat[,i] <-featureMatrix[,i]+sumvect
      }
    }





  if(length(label)==numSeqs){
    mismatchMat<-as.data.frame(mismatchMat)
    mismatchMat<-cbind(mismatchMat,label)
  }


  return(mismatchMat)
}


