#' K-Nearest Neighbor for Peptides (KNNPeptide)
#'
#'
#' This function needs an extra training data set and a label. We compute the similarity score of each input sequence with all sequences in the training data set.
#' We use the BLOSUM62 matrix to compute the similarity score. The label shows the class of each sequence in the training data set.
#' KNNPeptide finds the label of 1%...percent% of the most similar training sequence with the input sequence.
#' It reports the frequency of each class for each k% most similar sequences. The length of the feature vector will be percent*(number of classes).
#'
#' @note This function is usable for amino acid sequences with the same length in both training data set and the set of sequences.
#'
#' @references Chen, Zhen, et al. "iFeature: a python package and web server for features extraction and selection from protein and peptide sequences." Bioinformatics 34.14 (2018): 2499-2502.
#'
#' @param seqs is a fasta file with amino acids sequences. Each sequence starts
#' with a '>' character or it is a string vector such that each element is a peptide or protein sequence.
#'
#'
#' @param trainSeq is a fasta file with amino acids sequences. Each sequence starts
#' with a '>' character. Also it could be a string vector such that each element is a peptide sequence. Eaxh sequence in the training set
#' is associated with a label. The label is found in the parameret labeltr.
#'
#' @param labeltr This parameter is a vector whose length is equivalent to the number of sequences in the training set. It shows class of
#' each sequence in the trainig set.
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @param percent determines the threshold which is used to identify sequences (in the training set) which are similar to the input sequence.
#'
#' @return This function returns a feature matrix such that number of columns is number of classes multiplied by percent and number of rows is equal to the number of the sequences.
#'
#'
#' @export
#'
#' @examples
#'
#'
#' ptmSeqsADR<-system.file("extdata/",package="ftrCOOL")
#' ptmSeqsVect<-as.vector(read.csv(paste0(ptmSeqsADR,"/ptmVect101AA.csv"))[,2])
#'
#' posSeqs<-as.vector(read.csv(paste0(ptmSeqsADR,"/poSeqPTM101.csv"))[,2])
#' negSeqs<-as.vector(read.csv(paste0(ptmSeqsADR,"/negSeqPTM101.csv"))[,2])
#'
#' posSeqs<-posSeqs[1:10]
#' negSeqs<-negSeqs[1:10]
#'
#' trainSeq<-c(posSeqs,negSeqs)
#'
#' labelPos<-rep(1,length(posSeqs))
#' labelNeg<-rep(0,length(negSeqs))
#'
#' labeltr<-c(labelPos,labelNeg)
#'
#' KNNPeptide(seqs=ptmSeqsVect,trainSeq=trainSeq,percent=10,labeltr=labeltr)
#'



KNNPeptide <- function(seqs,trainSeq,percent=30,label=c(),labeltr=c())
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
    stop("ERROR, input sequence is not in a correct type. it should be a fasta file or a string vector.")
  }


  if(length(trainSeq)==1&&file.exists(trainSeq)){
    trainSeq<-fa.read(trainSeq,alphabet="aa")
    trseqs_Lab<-alphabetCheck(trainSeq,alphabet = "aa",label=labeltr)
    trainSeq<-trseqs_Lab[[1]]
    labeltr<-trseqs_Lab[[2]]
  }
  else if(is.vector(trainSeq)){

    trainSeq<-sapply(trainSeq,toupper)
    trseqs_Lab<-alphabetCheck(trainSeq,alphabet = "aa",label=labeltr)
    trainSeq<-trseqs_Lab[[1]]
    labeltr<-trseqs_Lab[[2]]
  }
  else {
    stop("ERROR, input sequence is not in a correct type. it should be a fasta file or a string vector.")
  }





  lenSeqs<-sapply(seqs,nchar)
  lenTrain<-sapply(seqs,nchar)
  len<-c(lenSeqs,lenTrain)
  len<-unique(len)

  if(length(len)>1){
    stop("ERROR all sequences should be in a same length both in train and sequences also train sequences and sample sequences has similar length")

  }
  numTrain<-length(trainSeq)
  if(length(labeltr)!=numTrain){
    stop("ERROR length label of trainSeq and number of sequences in train file should be the same")
  }

  #Blosum<-read.csv("R/BLOSUM62.csv")
  BlosumAD<-paste0(path.pack,"/BLOSUM62.csv")
  Blosum<-read.csv(BlosumAD)

  row.names(Blosum)<-Blosum[,1]
  Blosum<-Blosum[,-1]
  Blosum<-as.matrix(Blosum)
  Blosum<-type.convert(Blosum)
  numSeqs=length(seqs)

  charsTrain<-lapply(trainSeq, function(i) unlist(strsplit(i,"")))
  charsSeq<-lapply(seqs, function(i) unlist(strsplit(i,"")))
  simMatrix<-matrix(0,nrow = numSeqs,ncol = numTrain)

  mini=min(Blosum)
  maxi=max(Blosum)

  colnames(simMatrix)<-names(trainSeq)
  rownames(simMatrix)<-names(seqs)


  for(n in 1:numSeqs){
    charSeq<-charsSeq[[n]]
    for(m in 1:numTrain){
      charTrain<-charsTrain[[m]]
      blsum<-diag(Blosum[charSeq,charTrain])
      blsum[blsum<0]<-0


      simMatrix[n,m]<-sum((blsum-mini)/(maxi-mini))/len
    }
  }

  knum<-1:percent
  knum<-knum/100
  knum<-ceiling(knum*numTrain)

  #disMatrix<-t(apply(disMatrix, 1, sort))
  tabLabel<-table(labeltr)
  numClass<-length(tabLabel)
  featureMatrix<-vector()

  priorityLabel<-matrix(nrow = numSeqs,ncol = numTrain)
  #seqPeriority<-1:numTrain
  for(n in 1:numSeqs){
    temp<-t(rbind(simMatrix[n,],labeltr))
    temp<-temp[order(temp[,1]),]
    priorityLabel[n,]<-rev(temp[,2])

  }


  for(i in 1:length(knum)){

    classMatrix<-matrix(0,ncol = numClass,nrow = numSeqs)
    colnames(classMatrix)<-names(tabLabel)

    tempMatrix<-priorityLabel[,1:knum[i]]
    if(knum[i]>1&&numSeqs>1){
      for(j in 1:numSeqs){
        tvect<-tempMatrix[j,]
        tabtvect<-table(tvect)
        classMatrix[j,names(tabtvect)]<-tabtvect
      }
    }
    else {
      if(numSeqs==1){
        tableList<-table(tempMatrix)
        for(j in 1:numSeqs){
          classMatrix[j,names(tableList[[j]])]<-tableList[[j]]
        }
      }
      else{
        for(j in 1:numSeqs){
          classMatrix[j,as.character(tempMatrix)[j]]<-1
        }
      }

    }

    colnames(classMatrix)<-paste0("c",colnames(classMatrix),"_",i,"%")
    featureMatrix<-cbind(featureMatrix,classMatrix)
  }


  row.names(featureMatrix)<-names(seqs)
  return(featureMatrix)


}

