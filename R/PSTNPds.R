#' Position-Specific Trinucleotide Propensity based on double-strand (PSTNPds)
#'
#' This function works like \link{PSTNPss_DNA} except that it considers T as A and G as C. So it
#' converts Ts in the sequence to A and Gs to C. Then, it works with 2 alphabets A and C.
#' For more details refer to \link{PSTNPss_DNA}.
#'
#' @references Chen, Zhen, et al. "iLearn: an integrated platform and meta-learner for feature engineering, machine-learning analysis and modeling of DNA, RNA and protein sequence data." Briefings in bioinformatics 21.3 (2020): 1047-1057.
#'
#'
#' @note The length of the sequences in positive and negative data sets and the input sets
#' should be equal.
#'
#' @param seqs is a FASTA file containing nucleotide sequences. The sequences start
#' with '>'. Also, seqs could be a string vector. Each element of the vector is a nucleotide sequence.
#'
#' @param pos is a fasta file containing nucleotide sequences. Each sequence starts
#' with '>'. Also, the value of this parameter can be a string vector.
#' The sequences are positive sequences in the training model.
#'
#' @param neg is a fasta file containing nucleotide sequences. Each sequence starts
#' with '>'. Also, the value of this parameter can be a string vector.
#' The sequences are negative sequences in the training model.
#'
#' @param label is an optional parameter. It is a vector whose length is equal to the number of sequences.
#' It shows the class of each entry (i.e., sequence).
#'
#' @return It returns a feature matrix. The number of columns is equal to the length of sequences minus two
#' and the number of rows is equal to the number of sequences.
#'
#' @export
#'
#' @examples
#'
#' ptmSeqsADR<-system.file("extdata/",package="ftrCOOL")
#'
#' posSeqs<-fa.read(file=paste0(ptmSeqsADR,"/posData.txt"),alphabet="dna")
#' negSeqs<-fa.read(file=paste0(ptmSeqsADR,"/negData.txt"),alphabet="dna")
#' seqs<-fa.read(file=paste0(ptmSeqsADR,"/testData.txt"),alphabet="dna")
#'
#'
#' PSTNPds(seqs=seqs,pos=posSeqs[1],neg=negSeqs[1])
#'
#'
#'

PSTNPds<-function(seqs,pos,neg,label=c()){

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

  }else {
    stop("ERROR, input sequence is not in a correct type. It should be a FASTA file or a string vector.")
  }

  seqs<-gsub("T","A",seqs)
  seqs<-gsub("G","C",seqs)
  lens<-sapply(seqs,nchar)
  lenSeq<-unique(lens)
  if(length(lenSeq)>1){
    stop("Error sequences should be in the same length")
  }

  if(length(pos)==1&&file.exists(pos)){
    posSeqs<-fa.read(pos,alphabet="dna")
    posSeqs<-alphabetCheck(posSeqs,alphabet = "dna")

    posSeqs<-posSeqs[[1]]

  }
  else if(is.vector(pos)){
    posSeqs<-sapply(pos,toupper)

    posSeqs<-alphabetCheck(posSeqs,alphabet = "dna")

    posSeqs<-posSeqs[[1]]

  }else {
    stop("ERROR, positive sequences is not in a correct type. It should be a FASTA file or a string vector.")
  }

  posSeqs<-gsub("T","A",posSeqs)
  posSeqs<-gsub("G","C",posSeqs)

  lenPosSeqs<-sapply(posSeqs,nchar)
  lenPos<-unique(lenPosSeqs)
  if(length(lenPos)>1){
    stop("Error positive sequences should be in the same length")
  }

  if(lenPos!=lenSeq){
    stop("Posetive sequences and sample sequences should be in the same length")
  }

  if(length(neg)==1&&file.exists(neg)){
    negSeqs<-fa.read(neg,alphabet="dna")
    negSeqs<-alphabetCheck(negSeqs,alphabet = "dna")

    negSeqs<-negSeqs[[1]]

  }
  else if(is.vector(neg)){
    negSeqs<-sapply(neg,toupper)

    negSeqs<-alphabetCheck(negSeqs,alphabet = "dna")

    negSeqs<-negSeqs[[1]]

  }else {
    stop("ERROR, negative sequences is not in a correct type. It should be a FASTA file or a string vector.")
  }

  negSeqs<-gsub("T","A",negSeqs)
  negSeqs<-gsub("G","C",negSeqs)

  lenNegSeqs<-sapply(negSeqs,nchar)
  lenNeg<-unique(lenNegSeqs)

  if(length(lenNeg)>1){
    stop("Error negative sequences should be in the same length")
  }

  if(lenNeg!=lenSeq){
    stop("Error negative sequences and sample sequences should be in the same length")
  }



  tripletPos<-sapply(posSeqs, function(x) {temp<-unlist(strsplit(x,split = ""))
  len=length(temp)
  temp1<-temp[1:(len-2)]
  temp2<-temp[2:(len-1)]
  temp3<-temp[3:len]
  paste(temp1,temp2,temp3,sep = "")
  })

  tripletPos<-t(tripletPos)
  tabPos<-apply(tripletPos, 2, table)

  tripletNeg<-sapply(negSeqs, function(x) {temp<-unlist(strsplit(x,split = ""))
  len=length(temp)
  temp1<-temp[1:(len-2)]
  temp2<-temp[2:(len-1)]
  temp3<-temp[3:len]
  paste(temp1,temp2,temp3,sep = "")
  })

  tripletNeg<-t(tripletNeg)
  tabNeg<-apply(tripletNeg, 2, table)

  posMat<-matrix(0,ncol = (lenSeq-2),nrow = 8)
  negMat<-matrix(0,ncol = (lenSeq-2),nrow = 8)
  rNams<-nameKmer(k=3,type="num",2)
  rNams<-c("AAA","AAC","ACA","ACC","CAA","CAC","CCA","CCC")
  row.names(posMat)<-rNams
  row.names(negMat)<-rNams
  for(i in 1:(lenSeq-2)){
    if(length(posSeqs)>1){
      np=names(tabPos[[i]])
    } else{
      np=tripletPos[i]
    }

    posMat[np,i]=tabPos[[i]]

    if(length(posSeqs)>1){
      nn=names(tabNeg[[i]])
    } else{
      nn=tripletNeg[i]
    }
    negMat[nn,i]=tabNeg[[i]]
  }

  z<- posMat-negMat

  tripletSamples<-lapply(seqs, function(x) {temp<-unlist(strsplit(x,split = ""))
  len=length(temp)
  temp1<-temp[1:(len-2)]
  temp2<-temp[2:(len-1)]
  temp3<-temp[3:len]
  paste(temp1,temp2,temp3,sep = "")
  })

  tripletSamples<-tripletSamples
  outPutMat<-matrix(0,nrow = length(seqs),ncol=(lenSeq-2))
  row.names(outPutMat)<-names(seqs)
  for(n in 1:length(seqs)){
    outPutMat[n,]<-diag(z[tripletSamples[[n]],1:(lenSeq-2)])
  }



  return(outPutMat)

}


