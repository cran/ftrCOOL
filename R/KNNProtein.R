#' K-Nearest Neighbor for Protein (KNNProtein)
#'
#' This function is like \link{KNNPeptide} with the difference that similarity score is computed by Needleman-Wunsch algorithm.
#'
#' @references Chen, Zhen, et al. "iFeature: a python package and web server for features extraction and selection from protein and peptide sequences." Bioinformatics 34.14 (2018): 2499-2502.
#'
#' @param seqs is a fasta file with amino acids sequences. Each sequence starts
#' with a '>' character. Also it could be a string vector such that each element is a protein sequence.
#'
#' @param trainSeq is a fasta file with amino acids sequences. Each sequence starts
#' with a '>' character. Also it could be a string vector such that each element is a protein sequence. Eaxh sequence in the training set
#' is associated with a label. The label is found in the parameret labeltr.
#'
#' @param labeltr This parameter is a vector whose length is equivalent to the number of sequences in the training set. It shows class of
#' each sequence in the trainig set.
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
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
#' ptmSeqsADR<-system.file("extdata/",package="ftrCOOL")
#' ptmSeqsVect<-as.vector(read.csv(paste0(ptmSeqsADR,"/ptmVect101AA.csv"))[,2])
#' ptmSeqsVect<-ptmSeqsVect[1:2]
#'
#' posSeqs<-as.vector(read.csv(paste0(ptmSeqsADR,"/poSeqPTM101.csv"))[,2])
#' negSeqs<-as.vector(read.csv(paste0(ptmSeqsADR,"/negSeqPTM101.csv"))[,2])
#'
#' posSeqs<-posSeqs[1:3]
#' negSeqs<-negSeqs[1:3]
#'
#' trainSeq<-c(posSeqs,negSeqs)
#'
#' labelPos<-rep(1,length(posSeqs))
#' labelNeg<-rep(0,length(negSeqs))
#'
#' labeltr<-c(labelPos,labelNeg)
#'
#' KNNProtein(seqs=ptmSeqsVect,trainSeq=trainSeq,percent=10,labeltr=labeltr)



KNNProtein <- function(seqs,trainSeq,percent=30,labeltr=c(),label=c())
{

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

  numTrain<-length(trainSeq)
  if(length(labeltr)!=numTrain){
    stop("Error, number of labels is not compatible with number of sequences in the training set")  }

  numSeqs=length(seqs)
  # simMatrix<-matrix(0,nrow = numSeqs,ncol = numTrain)
  # colnames(simMatrix)<-names(trainSeq)
  # rownames(simMatrix)<-names(seqs)

  temp1<-rep(1:numSeqs,each=numTrain)
  temp2<-rep(1:numTrain,numSeqs)

  listSeq<-list()
  for(i in 1:length(temp1)){
    listSeq[[i]]<-c(seqs[temp1[i]],trainSeq[temp2[i]])
  }

  simlist<-lapply(listSeq, function(i){
    seq1=i[1]
    seq2=i[2]
    gap=-1
    mismatch=-1
    match=1
    # Stop conditions
    stopifnot(gap <= 0) # check if penalty negative
    stopifnot(mismatch <= 0)  # check if penalty negative
    stopifnot(match >= 0)  # check if score positive

    # Initialize col and rownames for matrices
    len1 = nchar(seq1); len2 = nchar(seq2) # Save number of chars in each sequence
    seq1 = unlist(strsplit(seq1, split = "")) # convert seq to character vector
    seq2 = unlist(strsplit(seq2, split = "")) # convert seq to character vector

    # Initialize matrix M (for scores)
    M = matrix(0, nrow = len1 + 1, ncol = len2 + 1) # Initialize matrix
    rownames(M) = c("-", seq1) # assign seq chars to matrix names
    colnames(M) = c("-", seq2) # assign seq chars to matrix names
    M[1, ] = cumsum(c(0, rep(gap, len2))) # Fill 1st row with gap penalites
    M[, 1] = cumsum(c(0, rep(gap, len1))) # Fill 1st col with gap penalites

    # Initialize matrix D (for directions)
    D = matrix(0, nrow = len1 + 1, ncol = len2 + 1) # Initialize matrix
    rownames(D) = c("-", seq1) # assign seq chars to matrix names
    colnames(D) = c("-", seq2) # assign seq chars to matrix names
    D[1, ] = rep("hor") # Fill 1st row with "hor" for horizontal moves
    D[, 1] = rep("ver") # Fill 1st col with "ver" for vertical moves
    type = c("dia", "hor", "ver") # Lookup vector

    # Compute scores and save moves
    for (i in 2:(len1 + 1)){# for every (initially zero) row
      for (j in 2:(len2 + 1)){# for every (initially zero) col
        hor = M[i, j - 1] + gap # horizontal move = gap for seq1
        ver = M[i - 1, j] + gap # vertical move = gap for seq2
        dia = ifelse(rownames(M)[i] == colnames(M)[j], # diagonal = ifelse(chars equal, match, mismatch)
                     M[i - 1, j - 1] + match,
                     M[i - 1, j - 1] + mismatch)
        M[i, j] = max(dia, hor, ver) # Save current (best) score in M
        D[i, j] = type[which.max(c(dia, hor, ver))] # Save direction of move in D
      }
    }


    return(M[nrow(M), ncol(M)])

  })

  simvect<-unlist(simlist)
  simMatrix<-matrix(simvect,byrow = TRUE,ncol=numTrain,nrow=numSeqs)
  colnames(simMatrix)<-names(trainSeq)
  rownames(simMatrix)<-names(seqs)



  knum<-1:percent
  knum<-knum/100
  knum<-ceiling(knum*numTrain)

  #simMatrix<-t(apply(simMatrix, 1, sort))
  tabLabel<-table(labeltr)
  numClass<-length(tabLabel)
  featureMatrix<-vector()


  priorityLabel<-matrix(nrow = numSeqs,ncol = numTrain)


  for(n in 1:numSeqs){
    temp<-t(rbind(simMatrix[n,],labeltr))
    temp<-temp[order(temp[,1]),]
    priorityLabel[n,]<-rev(temp[,2])
  }

  featureMatrix<-vector()


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
