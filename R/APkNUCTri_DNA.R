#' Amphiphilic Pseudo-k Nucleotide Composition-Tri(series) (APkNUCTri_DNA)
#'
#'
#' This function calculates the amphiphilic pseudo k nucleotide composition(Tri) (Series)
#' for each sequence.
#'
#'
#' @details This function computes the pseudo nucleotide composition for each physicochemical property of trinucleotides.
#' We have provided users with the ability to choose among the 12 properties in the tri-nucleotide index database.
#'
#' @param seqs is a FASTA file containing nucleotide sequences. The sequences start
#' with '>'. Also, seqs could be a string vector. Each element of the vector is a nucleotide sequence.
#'
#' @param selectedIdx is a vector of Ids or indices of the desired physicochemical properties of trinucleotides.
#' Users can choose the desired indices by their ids or their names in the TRI_DNA index file.
#' The default value of the parameter is a vector with ("Dnase I", "Bendability (DNAse)") ids.
#'
#' @param lambda is a tuning parameter. This integer value shows the maximum limit of spaces between trinucleotide pairs. The Number of spaces
#' changes from 1 to lambda.
#'
#' @param w (weight) is a tuning parameter. It changes in the range of 0 to 1. The default value is 0.05.
#'
#' @param l This parameter keeps the value of l in lmer composition. The lmers form the first 4^l of the APkNCTri descriptor.
#'
#' @param ORF (Open Reading Frame) is a logical parameter. If it is set to true, ORF region of each sequence is considered instead of the original sequence (i.e., 3-frame).
#'
#' @param reverseORF is a logical parameter. It is enabled only if ORF is true.
#' If reverseORF is true, ORF region will be searched in the sequence and also in the reverse complement of the sequence (i.e., 6-frame).
#'
#' @param threshold is a number between (0 , 1]. In selectedIdx, indices with a correlation
#' higher than the threshold will be deleted. The default value is 1.
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @return It is a feature matrix. The number of columns is 4^l+(number of the chosen indices*lambda) and the number of rows is equal to the number of sequences.
#'
#' @export
#'
#' @examples
#'
#' fileLNC<-system.file("extdata/Athaliana_LNCRNA.fa",package="ftrCOOL")
#' mat<-APkNUCTri_DNA(seqs=fileLNC,l=3,threshold=1)
#'

APkNUCTri_DNA<-function(seqs,selectedIdx=c("Dnase I", "Bendability (DNAse)"),lambda=3,w = 0.05,l=3,ORF=FALSE,reverseORF=TRUE,threshold=1,label=c())
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

  if(numSeqs==0){
    stop("There isn't any ORF region in these sequences")
  }




  #nucIdxAd<-paste0(path.pack,"/standardizeNUCidx3.csv")
  nucIdxAd<-paste0(path.pack,"/TRI_DNA.csv")
  nucIdx<-read.csv(nucIdxAd)

  row.names(nucIdx)<-nucIdx[,1]
  nucIdx<-nucIdx[selectedIdx,-1]
  nucIdx<-as.matrix(nucIdx)
  nucIdx<-type.convert(nucIdx)


  ###deleteing high correlate features
  if(threshold!=1){

    nucIdx<-t(nucIdx)
    corr<-cor(nucIdx)
    corr2<-corr^2
    tmp<-corr2
    tmp[upper.tri(tmp)]<-0
    for(i in 1:ncol(tmp)){
      tmp[i,i]=0
    }
    nucIdx<- nucIdx[,!apply(tmp,2,function(x) any(x > threshold))]
    nucIdx<- t(nucIdx)
  }
  ##normalize beween [0,1]
  minFea<-apply(nucIdx, 1, min)
  maxFea<-apply(nucIdx, 1, max)
  nucIdx<-(nucIdx-minFea)/(maxFea-minFea)

  numFea<-nrow(nucIdx)
  start<-4^l+1
  end<-4^l+(lambda*numFea)
  featureMatrix<-matrix(0,nrow = numSeqs,ncol = end)
  N<-sapply(seqs,nchar)

  TauMat<-matrix(0,nrow = numSeqs,ncol = (lambda*numFea))
  for(n in 1:numSeqs){
    seq<-seqs[n]

    if (lambda>N[n] || lambda<=0){
      stop("ERROR: lambda should be between [1,N]. N is the minimum of sequence lengths")
    }


    Tau<-vector(mode = "numeric")
    chars<-unlist(strsplit(seq,NULL))
    temp1<-chars[1:(N[n]-2)]
    temp2<-chars[2:(N[n]-1)]
    temp3<-chars[3:N[n]]
    Trimers<-paste0(temp1,temp2,temp3)
    lenTrimer=N[n]-2

    for(k in 1:lambda){
      vecti=1:(lenTrimer-k)
      vectj=vecti+k
      tempMat<-matrix(0,nrow = numFea,ncol = (lenTrimer-k))
      for(m in 1:numFea){
        tempMat[m,]=(nucIdx[m,Trimers[vecti]]*nucIdx[m,Trimers[vectj]])
      }
      Hi<-apply(tempMat, 1, sum)

      Tau<-c(Tau,Hi)
    }
    TauMat[n,]<-Tau

  }

  sum_Tau<-apply(TauMat, 1, sum)
  NUCmat<-kNUComposition_DNA(seqs,rng=l,normalized = FALSE,upto = FALSE)
  index=1
  for(index in 1:numSeqs){
    featureMatrix[index,1:(4^l)]<-NUCmat[index,]/(N[index]+(w*sum_Tau[index]))
    featureMatrix[index,start:end]<- w*TauMat[index,]/(N[index]+w*(sum_Tau[index]))


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
