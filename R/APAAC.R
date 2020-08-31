#' Amphiphilic Pseudo-Amino Acid Composition(APAAC) (series)
#'
#' This function calculates the amphiphilic pseudo amino acid composition (Series)
#' for each sequence.
#'
#' @details This function computes the pseudo amino acid composition for each physicochemical property.
#' We have provided users with the ability to choose among different properties (i.e., not confined to hydrophobicity or hydrophilicity).
#'
#'
#' @param seqs is a FASTA file with amino acid sequences. Each sequence starts
#' with a '>' character. Also, seqs could be a string vector. Each element of the vector is a peptide/protein sequence.
#'
#' @param aaIDX is a vector of Ids or indexes of the user-selected physicochemical properties in the aaIndex2 database.
#' The default values of the vector are the hydrophobicity ids and hydrophilicity ids
#' in the amino acid index file.
#'
#' @param lambda is a tuning parameter. Its value indicates the maximum number of spaces between di-nucleotide pairs. The number
#' changes from 1 to lambda.
#'
#' @param w (weight) is a tuning parameter. It changes in from 0 to 1. The default value is 0.5.
#'
#' @param l This parameter keeps the value of l in lmer composition. The lmers form the first 20^l elements of the APAAC descriptor.
#'
#' @param threshold is a number between (0 , 1]. In aaIDX, indices with a correlation
#' higher than the threshold will be deleted. The default value is 1.
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#'
#' @return A feature matrix such that the number of columns is 20^l+(number of chosen aaIndex*lambda) and the number of rows equals the number of sequences.
#'
#'
#' @export
#'
#' @examples
#'
#' filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
#' mat<-APAAC(seqs=filePrs,l=2,lambda=3,threshold=1)
#'

APAAC<-function(seqs,aaIDX=c("ARGP820101","HOPT810101"),lambda=3,w = 0.5,l=1,threshold=1,label=c())
{

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




  aaIdxAd<-paste0(path.pack,"/standardizedAAindex_555.csv")
  aaIdx<-read.csv(aaIdxAd)

  row.names(aaIdx)<-aaIdx[,1]
  aaIdx<-aaIdx[aaIDX,-1]
  aaIdx<-as.matrix(aaIdx)
  aaIdx<-type.convert(aaIdx)


  ###deleteing high correlated features
  if(threshold!=1){

    aaIdx<-t(aaIdx)
    corr<-cor(aaIdx)
    corr2<-corr^2
    tmp<-corr2
    tmp[upper.tri(tmp)]<-0
    for(i in 1:ncol(tmp)){
      tmp[i,i]=0
    }
    TFidx<-apply(tmp,2,function(x) any(x > threshold))
    DlIDX<-which(TFidx==TRUE)
    RMIDX<-which(TFidx==FALSE)


    if(length(DlIDX)!=0){
      delPropName<-names(DlIDX)
      delPropName<-toString(delPropName)
      warMessage<-paste("The sequences (",delPropName,") were deleted. They had a correlation with other properties more than the threshold")
      message(warMessage)
    }

    aaIdx<- aaIdx[,RMIDX]
    aaIdx<- t(aaIdx)
  }
  ##normalize beween [0,1]
  minFea<-apply(aaIdx, 1, min)
  maxFea<-apply(aaIdx, 1, max)
  aaIdx<-(aaIdx-minFea)/(maxFea-minFea)

  numFea<-nrow(aaIdx)
  start<-20^l+1
  end<-20^l+(lambda*numFea)
  featureMatrix<-matrix(0,nrow = numSeqs,ncol = end)
  for(n in 1:numSeqs){
    seq<-seqs[n]
    N<-nchar(seq)
    if (lambda>N || lambda<=0){
      stop("ERROR: lambda should be between [1,N]. N is the minimum of sequence lengths")
    }


    Tau<-vector(mode = "numeric")
    chars<-unlist(strsplit(seq,NULL))

    for(k in 1:lambda){
      vecti=1:(N-k)
      vectj=vecti+k
      tempMat<-matrix(0,nrow = numFea,ncol = (N-k))
      for(m in 1:numFea){
        tempMat[m,]=(aaIdx[m,chars[vecti]]*aaIdx[m,chars[vectj]])
      }
      Hi<-apply(tempMat, 1, sum)


      Tau<-c(Tau,Hi)
    }
    sum_Tau<-sum(Tau)
    AAC<-kAAComposition(seqs=seq,rng=l,normalized = FALSE,upto = FALSE)
    featureVector <- vector(mode = "numeric", length = length(Tau))
    featureVector[1:length(AAC)]<-(AAC)/(N+w*(sum_Tau))

    featureVector[start:end]<- (w*Tau)/(N+w*(sum_Tau))


    featureMatrix[n,]=featureVector


  }
  namAAC<-nameKmer(k=l,type = "aa")
  temp1=rep(1:lambda,each=numFea)
  temp2=rep(1:numFea,lambda)
  temp2=paste0("Fea",temp2)
  temp1=paste0("lambda",temp1)
  temp3=paste(temp1,temp2)
  colnames(featureMatrix)=c(namAAC,temp3)
  if(length(label)==numSeqs){
    featureMatrix<-as.data.frame(featureMatrix)
    featureMatrix<-cbind(featureMatrix,label)
  }
  row.names(featureMatrix)<-names(seqs)
  return(featureMatrix)
}
