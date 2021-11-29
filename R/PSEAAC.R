#' Pseudo-Amino Acid Composition (Parallel) (PSEAAC)
#'
#' This function calculates the pseudo amino acid composition (parallel)
#' for each sequence.
#'
#'
#' @param seqs is a FASTA file with amino acid sequences. Each sequence starts
#' with a '>' character. Also, seqs could be a string vector. Each element of the vector is a peptide/protein sequence.
#'
#' @param aaIDX is a vector of Ids or indexes of the user-selected physicochemical properties in the aaIndex2 database.
#' The default values of the vector are the hydrophobicity ids and hydrophilicity ids and Mass of residual
#' in the amino acid index file.
#'
#'
#' @param lambda is a tuning parameter. Its value indicates the maximum number of spaces between amino acid pairs. The number
#' changes from 1 to lambda.
#'
#' @param w (weight) is a tuning parameter. It changes in from 0 to 1. The default value is 0.05.
#'
#' @param l This parameter keeps the value of l in lmer composition. The lmers form the first 20^l elements of the APAAC descriptor.
#'
#' @param threshold is a number between (0 , 1]. It deletes aaIndexes which have a correlation
#' bigger than the threshold. The default value is 1.
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#'
#' @return A feature matrix such that the number of columns is 20^l+(lambda) and the number of rows is equal to the number of sequences.
#'
#' @export
#'
#' @examples
#' filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
#' mat<-PSEAAC(seqs=filePrs,l=2)


PSEAAC<-function(seqs,aaIDX=c("ARGP820101","HOPT810101","Mass"),lambda=30,w = 0.05,l=1,threshold=1,label=c()){

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

  numSeqs<-length(seqs)

   aaIdxAD<-paste0(path.pack,"/standardizedAAindex_555.csv")
   aaIdx<-read.csv(aaIdxAD)
   row.names(aaIdx)<-aaIdx[,1]
   aaIdx<-aaIdx[aaIDX,-1]
   aaIdx<-as.matrix(aaIdx)
   aaIdx<-type.convert(aaIdx)


  ###deleteing high correlate features
  if(threshold!=1){

    aaIdx<-t(aaIdx)
    corr<-cor(aaIdx)
    corr2<-corr^2
    tmp<-corr2
    tmp[upper.tri(tmp)]<-0
    for(i in 1:ncol(tmp)){
      tmp[i,i]=0
    }
    aaIdx<- aaIdx[,!apply(tmp,2,function(x) any(x > threshold))]
    aaIdx<- t(aaIdx)
  }##normalize beween [0,1]
  minFea<-apply(aaIdx, 1, min)
  maxFea<-apply(aaIdx, 1, max)
  aaIdx<-(aaIdx-minFea)/(maxFea-minFea)

  numFea<-nrow(aaIdx)

  start<-20^l+1
  end<-20^l+lambda

  featureMatrix<-matrix(0,nrow = numSeqs,ncol = ((20^l)+lambda))
  for(n in 1:numSeqs){
    seq<-seqs[n]
    N<-nchar(seq)
    if (lambda>N || lambda<=0){
      stop("Error: lambda should be between [1,N]. N is the minimum of sequence lengths")
    }


    small_theta<-vector(mode = "numeric", length = lambda)
    chars<-unlist(strsplit(seq,NULL))

    for(k in 1:lambda){
      vecti=1:(N-k)
      vectj=vecti+k
      tempMat<-matrix(0,nrow = numFea,ncol = (N-k))
      for(m in 1:numFea){
        tempMat[m,]=(aaIdx[m,chars[vecti]]-aaIdx[m,chars[vectj]])^2
      }
      bigThetaVect<-apply(tempMat, 2, sum)
      sumbigThetas<-sum(bigThetaVect)
      small_theta[k]<-(1/(N-k))*sumbigThetas
    }

    sum_small_th<-sum(small_theta)
    AAC<-kAAComposition(seq,rng=l,normalized = FALSE,upto = FALSE)
    featureVector <- vector(mode = "numeric", length = (length(AAC)+lambda))
    featureVector[1:length(AAC)]<-(AAC)/(N+w*(sum_small_th))

    featureVector[start:end]<- (w*small_theta)/(N+w*(sum_small_th))


    featureMatrix[n,]=featureVector


  }
  namAAC<-nameKmer(k=l,type = "aa")
  temp=1:lambda
  colnam=c(namAAC,paste("lambda",temp,sep=""))
  colnames(featureMatrix)=colnam
  if(length(label)==numSeqs){
    featureMatrix<-as.data.frame(featureMatrix)
    featureMatrix<-cbind(featureMatrix,label)
  }
  row.names(featureMatrix)<-names(seqs)

  return(featureMatrix)

}



