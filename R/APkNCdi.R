#' Amphiphilic Pseudo-k Nucleotide Composition-di(series)
#'
#'
#' This function calculates the amphiphilic pseudo k nucleotide composition(Di) (Series)
#' for each sequence.
#'
#'
#' @details This function computes the pseudo nucleotide composition for each physicochemical property of dinucleotides.
#' We have provided users with the ability to choose among the 125 properties in the di-nucleotide index database.
#'
#' @param seqs is a FASTA file containing nucleotide sequences. The sequences start
#' with '>'. Also, seqs could be a string vector. Each element of the vector is a nucleotide sequence.
#'
#' @param selectedNucIdx is a vector of Ids or indices of the desired physicochemical properties of dinucleotides.
#' Users can choose the desired indices by their ids or their names in the di-nucleotide index file.
#'
#' @param lambda is a tuning parameter. This integer value shows the maximum limit of spaces between dinucleotide pairs. The Number of spaces
#' changes from 1 to lambda.
#'
#' @param w (weight) is a tuning parameter. It changes in the range of 0 to 1. The default value is 0.5.
#'
#' @param l This parameter keeps the value of l in lmer composition. The lmers form the first 4^l elements of the APkNCdi descriptor.
#'
#' @param ORF (Open Reading Frame) is a logical parameter. If it is set to true, ORF region of each sequence is considered instead of the original sequence (i.e., 3-frame).
#'
#' @param reverseORF is a logical parameter. It is enabled only if ORF is true.
#' If reverseORF is true, ORF region will be searched in the sequence and also in the reverse complement of the sequence (i.e., 6-frame).
#'
#' @param threshold is a number between (0 to 1]. In selectedNucIdx, indices with a correlation
#' higher than the threshold will be deleted. The default value is 1.
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @return It is a feature matrix. The number of columns is 4^l+(number of the chosen indices*lambda) and the number of rows is equal to the number of sequences.
#'
#'
#' @export
#'
#' @examples
#'
#' fileLNC<-system.file("extdata/Athaliana_LNCRNA.fa",package="ftrCOOL")
#' mat<-APkNUCdi(seqs=fileLNC,selectedNucIdx=1:125,ORF=TRUE,threshold=0.8)
#'


APkNUCdi<-function(seqs,selectedNucIdx,lambda=3,w = 0.05,l=2,ORF=FALSE,reverseORF=TRUE,threshold=1,label=c())
{

  path.pack=system.file("extdata",package="ftrCOOL")

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
  flag=0
  if(ORF==TRUE){
    if(length(label)==length(seqs)){
      names(label)=names(seqs)
      flag=1
    }
    seqs=maxORF(seqs,reverse=reverseORF)
    if(flag==1)
      label=label[names(seqs)]
  }
  numSeqs<-length(seqs)



  nucIdxAd<-paste0(path.pack,"/standardizeNUCidx2.csv")
  nucIDX<-read.csv(nucIdxAd)

  row.names(nucIDX)<-nucIDX[,1]
  nucIDX<-nucIDX[selectedNucIdx,-1]
  nucIDX<-as.matrix(nucIDX)
  nucIDX<-type.convert(nucIDX)


  ###deleteing high correlate features
  if(threshold!=1){

    nucIDX<-t(nucIDX)
    corr<-cor(nucIDX)
    corr2<-corr^2
    tmp<-corr2
    tmp[upper.tri(tmp)]<-0
    for(i in 1:ncol(tmp)){
      tmp[i,i]=0
    }
    nucIDX<- nucIDX[,!apply(tmp,2,function(x) any(x > threshold))]
    nucIDX<- t(nucIDX)
  }
  ##normalize beween [0,1]
  minFea<-apply(nucIDX, 1, min)
  maxFea<-apply(nucIDX, 1, max)
  nucIDX<-(nucIDX-minFea)/(maxFea-minFea)

  numFea<-nrow(nucIDX)
  start<-4^l+1
  end<-4^l+(lambda*numFea)
  featureMatrix<-matrix(0,nrow = numSeqs,ncol = end)
  for(n in 1:numSeqs){
    seq<-seqs[n]
    N<-nchar(seq)
    if (lambda>N || lambda<=0){
      stop("ERROR: lambda should be between [1,N]. N is the minimum of sequence lengths")
    }


    Tau<-vector(mode = "numeric")
    chars<-unlist(strsplit(seq,NULL))
    temp1<-chars[1:(N-1)]
    temp2<-chars[2:N]
    dimers<-paste0(temp1,temp2)
    lenDimer=N-1

    for(k in 1:lambda){
      vecti=1:(lenDimer-k)
      vectj=vecti+k
      tempMat<-matrix(0,nrow = numFea,ncol = (lenDimer-k))
      for(m in 1:numFea){
        tempMat[m,]=(nucIDX[m,dimers[vecti]]*nucIDX[m,dimers[vectj]])
      }
      Hi<-apply(tempMat, 1, sum)


      Tau<-c(Tau,Hi)
    }
    sum_Tau<-sum(Tau)
    NUC<-kNUComposition(seq,rng=l,normalized = FALSE,upto = FALSE)
    featureVector <- vector(mode = "numeric", length = length(Tau))
    featureVector[1:length(NUC)]<-(NUC)/(N+w*(sum_Tau))

    featureVector[start:end]<- (w*Tau)/(N+w*(sum_Tau))


    featureMatrix[n,]=featureVector


  }
  namNUC<-nameKmer(k=l,type = "dna")
  temp1=rep(1:lambda,each=numFea)
  temp2=rep(1:numFea,lambda)
  temp2=paste0("Fea",temp2)
  temp1=paste0("lambda",temp1)
  temp3=paste(temp1,temp2)
  colnames(featureMatrix)=c(namNUC,temp3)
  if(length(label)==numSeqs){
    featureMatrix<-as.data.frame(featureMatrix)
    featureMatrix<-cbind(featureMatrix,label)
  }
  row.names(featureMatrix)<-names(seqs)
  return(featureMatrix)
}
