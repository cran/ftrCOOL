#' Novel Pseudo k Nucleotide Composition (series)
#'
#' This function replaces nucleotides with a four-length vector.
#' The first three elements represent the nucleotides and
#' the forth holds the frequency of the nucleotide from the beginning of the sequence until the position of the nucleotide in the sequence.
#' 'A' will be replaced with c(1, 1, 0, freq), 'C' with c(0, 1, 1, freq),'G' with c(1, 0, 1, freq), and 'T' with c(0, 0, 0, freq).
#'
#' @references Feng, Pengmian, et al. "iDNA6mA-PseKNC: Identifying DNA N6-methyladenosine sites by incorporating nucleotide physicochemical properties into PseKNC." Genomics 111.1 (2019): 96-102.
#'
#' @param seqs is a FASTA file containing nucleotide sequences. The sequences start
#' with '>'. Also, seqs could be a string vector. Each element of the vector is a nucleotide sequence.
#'
#'
#' @param outFormat (output format) can take two values: 'mat'(matrix) and 'txt'. The default value is 'mat'.
#'
#' @param outputFileDist shows the path and name of the 'txt' output file.
#'
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @return A feature matrix. The number of rows is equal to the number of sequences.
#'
#'
#' @export
#'
#' @examples
#'
#' dir = tempdir()
#' LNCSeqsADR<-system.file("extdata/",package="ftrCOOL")
#' LNC50Nuc<-as.vector(read.csv(paste0(LNCSeqsADR,"/LNC50Nuc.csv"))[,2])
#' mat<-novel_PseKNC(seqs = LNC50Nuc,outFormat="mat")
#'
#' ad<-paste0(dir,"/ENUCcompos.txt")
#' fileLNC<-system.file("extdata/Athaliana_LNCRNA.fa",package="ftrCOOL")
#' novel_PseKNC(seqs = fileLNC,outFormat="txt",outputFileDist=ad)

novel_PseKNC<-function(seqs,outFormat="mat",outputFileDist="",label=c()){


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



  nucs<-list("A"=c(1,1,0),"C"=c(0,1,1),"G"=c(1,0,1),"T"=c(0,0,0))
  numSeqs<-length(seqs)

  if(outFormat=="mat"){

    if(length(unique(lenSeqs))>1){
      stop("ERROR: All sequences should have the same length in 'mat' mode. For sequences with different lengths, please use 'txt' for outFormat parameter")
    }
  featureMatrix<-matrix(0,nrow = numSeqs,ncol = (4*lenSeqs[1]))
  for(n in 1:numSeqs){
    seq<-seqs[n]
    N<-nchar(seq)
    freq<-vector()
    featureVect<-vector()
    charList<-unlist(strsplit(seq,split = ""))
    for(i in 1:N){
      a=freq[charList[i]]
      if(is.na(a)){
        freq<-c(freq,1)
        names(freq)[length(freq)]<-charList[i]
        featureVect<-c(featureVect,c(nucs[[charList[i]]],freq[charList[i]]/i))

      } else {
        freq[charList[i]]<-freq[charList[i]]+1
        featureVect<-c(featureVect,c(nucs[[charList[i]]],freq[charList[i]]/i))

      }

    }
    featureMatrix[n,]<-featureVect
  }
  temp1=rep(1:lenSeqs[1],each=4)
  temp1<-paste0("pos",temp1)
  temp2<-rep(c("Ring","Fun_A","H_bound_St","Freq"),lenSeqs[1])
  colnames(featureMatrix)<-paste0(temp1,temp2)
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
      N<-nchar(seq)
      freq<-vector()
      featureVect<-vector()
      charList<-unlist(strsplit(seq,split = ""))
      for(i in 1:N){
        a=freq[charList[i]]
        if(is.na(a)){
          freq<-c(freq,1)
          names(freq)[length(freq)]<-charList[i]
          featureVect<-c(featureVect,c(nucs[[charList[i]]],freq[charList[i]]/i))

        } else {
          freq[charList[i]]<-freq[charList[i]]+1
          featureVect<-c(featureVect,c(nucs[[charList[i]]],freq[charList[i]]/i))

        }

      }
      temp<-c(nameSeq[n],featureVect)
      temp<-paste(temp,collapse = "\t")
      write(temp,outputFileDist,append = TRUE)
    }
  }
  else {
    stop("ERROR: outFormat should be 'mat' or 'txt' ")
  }

}
