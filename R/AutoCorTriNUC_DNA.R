#' Tri Nucleotide Autocorrelation-Autocovariance (AutoCorTriNUC_DNA)
#'
#'
#' It creates the feature matrix for each function in autocorelation
#' (i.e., Moran, Greay, NormalizeMBorto) or autocovariance (i.e., AC, CC,ACC). The user can
#' select any combination of the functions too. In this case, the final matrix will contain
#' features of each selected function.
#'
#' @details For CC and AAC autocovriance functions, which consider the covariance of the
#' two physicochemical properties, we have provided users with the ability to categorize
#' their selected properties in a list.
#' The binary combination of each group will be taken into account.
#' Note: If all the features are in a group or selectedAAidx parameter is a vector,
#' the binary combination will be calculated for all the physicochemical properties.
#'
#'
#' @param seqs is a FASTA file containing nucleotide sequences. The sequences start
#' with '>'. Also, seqs could be a string vector. Each element of the vector is a nucleotide sequence.
#'
#'
#' @param selectedNucIdx function takes as input the physicochemical properties. Users select the properties by their ids
#' or indices in the TRI_DNA file. This parameter could be a vector or a list of trinucleotide indices.
#' The default value of this parameter is a vector with  ("Dnase I", "Bendability (DNAse)") ids.
#'
#'
#' @param maxlag This parameter shows the maximum gap between two tri-nucleotide pairs. The gaps change from 1 to maxlag (the maximum lag).
#'
#' @param type could be 'Moran', 'Greay', 'NormalizeMBorto', 'AC', 'CC', or 'ACC'. Also, it could be any combination of them.
#'
#' @param threshold is a number between (0 to 1]. In selectedNucIdx, indices with a correlation higher than the threshold will be deleted.The default value is 1.
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#'
#' @return This function returns a feature matrix. The number of columns in the matrix changes depending on the chosen autocorrelation or autocovariance types and nlag parameter.
#' The output is a matrix. The number of rows shows the number of sequences.
#'
#'
#' @export
#'
#' @examples
#'
#' fileLNC<-system.file("extdata/Athaliana_LNCRNA.fa",package="ftrCOOL")
#' mat1<-AutoCorTriNUC_DNA(seqs=fileLNC,selectedNucIdx=c(1:7),maxlag=20,type=c("Moran","Geary"))
#'
#' mat2<-AutoCorTriNUC_DNA(seqs=fileLNC,selectedNucIdx=list(c(1,3),6:10,c(2:7)),
#' maxlag=15,type=c("AC","CC"))
#'

AutoCorTriNUC_DNA<-function(seqs,selectedNucIdx=c("Dnase I", "Bendability (DNAse)"),maxlag=3,threshold=1,type=c("Moran","Geary","NormalizeMBorto","AC","CC","ACC"),label=c()){


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


  numSeqs<-length(seqs)
  lenSeqs<-sapply(seqs, nchar)
  lenSeqs<-lenSeqs-2
  Trimer<-lapply(seqs, function(i) {
    chars<-unlist(strsplit(i,split = ""))
    N=nchar(i)
    temp1<-chars[1:(N-2)]
    temp2<-chars[2:(N-1)]
    temp3<-chars[3:N]
    trimers<-paste0(temp1,temp2,temp3)
    return(trimers)
  })



  #aaIdxAd<-paste0(path.pack,"/standardizeNUCidx3.csv")
  aaIdxAd<-paste0(path.pack,"/TRI_DNA.csv")
  aaIdx<-read.csv(aaIdxAd)


  rNumAAIdx<-unlist(selectedNucIdx)
  unqrNumAAIdx<-unique(rNumAAIdx)
  row.names(aaIdx)<-aaIdx[,1]
  aaIdx<-as.matrix(aaIdx)

  selectedNucIdxNames=list()
  if(is.list(selectedNucIdx))
  {
    for(i in 1:length(selectedNucIdx)){
      if(is.numeric(selectedNucIdx[[i]])){
        names(selectedNucIdx[[i]])<-row.names(aaIdx)[selectedNucIdx[[i]]]
      }
      else if(is.character(selectedNucIdx[[i]])){
        names(selectedNucIdx[[i]])<-selectedNucIdx[[i]]
      } else{
        stop('ERROR: selecteNucIdx should be Ids or indices of the desired properties in the tri-nucleotide index file')
      }
    }

    selectedNucIdxNames<-lapply(selectedNucIdx, function(x){
      names(x)
    })
  }
  else if(is.vector(selectedNucIdx)){
    if(is.numeric(selectedNucIdx)){
      names(selectedNucIdx)<-row.names(aaIdx)[selectedNucIdx]
    } else if(is.character(selectedNucIdx)){
      names(selectedNucIdx)<-selectedNucIdx
    } else{
      stop('ERROR: selecteNucIdx should be Ids or indices of the desired properties in the tri-nucleotide index file')
    }
    selectedNucIdxNames=c()
    selectedNucIdxNames<-names(selectedNucIdx)
  }
  else{
    stop("ERROR: selectedAAidx should be a vector or a list")
  }
  aaIdx<-aaIdx[unqrNumAAIdx,-1]
  aaIdx<-type.convert(aaIdx)



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
      delPropNameSTR<-toString(delPropName)
      warMessage<-paste("The properties (",delPropNameSTR,") were deleted. They had a correlation with other properties more than the threshold")
      message(warMessage)
      selectedNucIdxNames<-lapply(selectedNucIdxNames, function(x){
        deletedIdx<-which(x %in% delPropName)
        if(length(deletedIdx)!=0){
          x=x[-deletedIdx]
        }
        return(x)
      })
    }


    aaIdx<- aaIdx[,RMIDX]
    aaIdx<- t(aaIdx)

  }

  numAAidx<-nrow(aaIdx)



  ## average all tri-nucleotides in seq for each property
  ##is a vector with length property number
  Pprim <- matrix(0, ncol = numAAidx,nrow = numSeqs)


  # Is a matrix each row belong to a sequence each col belong to each property
  for( j in 1:numSeqs){
    for(i in 1:numAAidx)
    {
      Pprim[j,i]<-sum(aaIdx[i,Trimer[[j]]])/lenSeqs[j]
    }
  }

  Pprim<-t(Pprim)
  row.names(Pprim)<-row.names(aaIdx)
  totalfeatureMatrix=vector()
  if("Moran" %in% type){
    featureMatrix<-Moran_CorTri(maxLag = maxlag,aaIdx = aaIdx,pPrim = Pprim,numSeqs=numSeqs,lenSeqs=lenSeqs,Trimer=Trimer)
    colnames(featureMatrix)<-paste("Mo",colnames(featureMatrix))
    totalfeatureMatrix<-cbind(totalfeatureMatrix,featureMatrix)
  }
  if("Geary" %in% type){
    featureMatrix<-Geary_CorTri(maxLag = maxlag,aaIdx = aaIdx,pPrim = Pprim,numSeqs=numSeqs,lenSeqs=lenSeqs,Trimer=Trimer)
    colnames(featureMatrix)<-paste("Ge",colnames(featureMatrix))
    totalfeatureMatrix<-cbind(totalfeatureMatrix,featureMatrix)
  }
  if("NormalizeMBorto" %in% type){
    featureMatrix<-NormalizeMBorto_CorTri(maxLag = maxlag,aaIdx = aaIdx,pPrim = Pprim,numSeqs=numSeqs,lenSeqs=lenSeqs,Trimer=Trimer)
    colnames(featureMatrix)<-paste("No",colnames(featureMatrix))
    totalfeatureMatrix<-cbind(totalfeatureMatrix,featureMatrix)
  }
  if("AC" %in% type){
    featureMatrix<-AC_CorTri(maxLag = maxlag,aaIdx = aaIdx,pPrim = Pprim,numSeqs=numSeqs,lenSeqs=lenSeqs,Trimer=Trimer)
    colnames(featureMatrix)<-paste("AC",colnames(featureMatrix))
    totalfeatureMatrix<-cbind(totalfeatureMatrix,featureMatrix)
  }
  if("CC" %in% type){
    featureMatrix<-CC_CorTri(selectedProp = selectedNucIdxNames ,maxLag = maxlag,aaIdx = aaIdx,pPrim = Pprim,numSeqs=numSeqs,lenSeqs=lenSeqs,Trimer=Trimer)
    colnames(featureMatrix)<-paste("CC",colnames(featureMatrix))
    totalfeatureMatrix<-cbind(totalfeatureMatrix,featureMatrix)
  }
  if("ACC" %in% type){
    featureMatrix<-ACC_CorTri(selectedProp = selectedNucIdxNames ,maxLag = maxlag,aaIdx = aaIdx,pPrim = Pprim,numSeqs=numSeqs,lenSeqs=lenSeqs,Trimer=Trimer)
    colnames(featureMatrix)<-paste("ACC",colnames(featureMatrix))
    totalfeatureMatrix<-cbind(totalfeatureMatrix,featureMatrix)
  }
  if(length(label)==numSeqs){
    featureMatrix<-as.data.frame(featureMatrix)
    featureMatrix<-cbind(featureMatrix,label)
  }
  row.names(featureMatrix)<-names(seqs)


  return(totalfeatureMatrix)

}

Moran_CorTri<-function(maxLag=30,aaIdx,pPrim,numSeqs,lenSeqs,Trimer){

  totalMatrix<-matrix(0,nrow = numSeqs,ncol = (nrow(aaIdx)*maxLag))
  tempName<-rep(1:maxLag,nrow(aaIdx))
  tempName<-paste0("d=",tempName)
  temp2Name<-rep(row.names(aaIdx),each=maxLag)
  colnames(totalMatrix)<-paste(tempName,temp2Name)
  for(n in 1:numSeqs){
    triNuc<-Trimer[[n]]
    N<-lenSeqs[n]
    featureMartix<-matrix(ncol = nrow(aaIdx), nrow = maxLag)
    tempMat<-(aaIdx[,triNuc]-pPrim[,n])^2
    sigmaRow<-apply(tempMat,1,sum)
    denominator<-(1/lenSeqs[n])*sigmaRow
    featVect=vector()
    for(d in 1:maxLag){
      vecti<-1:(N-d)
      vectj<-vecti+d
      tempMat<-(aaIdx[,triNuc[vecti]]-pPrim[,n])*(aaIdx[,triNuc[vectj]]-pPrim[,n])
      sigma<-apply(tempMat, 1, sum)
      numerator<-sigma*(1/(N-d))/denominator
      featVect<-rbind(featVect,numerator)
    }
    featVect<-as.vector(featVect)

    totalMatrix[n,]<-featVect

  }
  return(totalMatrix)
}

Geary_CorTri<-function(maxLag=30,aaIdx,pPrim,numSeqs,lenSeqs,Trimer){
  totalMatrix<-matrix(0,nrow = numSeqs,ncol = (nrow(aaIdx)*maxLag))
  tempName<-rep(1:maxLag,nrow(aaIdx))
  tempName<-paste0("d=",tempName)
  temp2Name<-rep(row.names(aaIdx),each=maxLag)
  colnames(totalMatrix)<-paste(tempName,temp2Name)
  for(n in 1:numSeqs){
    triNuc<-Trimer[[n]]
    N<-lenSeqs[n]
    featureMartix<-matrix(ncol = nrow(aaIdx), nrow = maxLag)
    tempMat<-(aaIdx[,triNuc]-pPrim[,n])^2
    sigmaRow<-apply(tempMat,1,sum)
    denominator<-(1/lenSeqs[n])*sigmaRow
    featVect=vector()
    for(d in 1:maxLag){
      vecti<-1:(N-d)
      vectj<-vecti+d
      tempMat<-(aaIdx[,triNuc[vecti]]-aaIdx[,triNuc[vectj]])^2
      sigma<-apply(tempMat, 1, sum)
      numerator<-sigma*(1/(N-d))/denominator
      featVect<-rbind(featVect,numerator)
    }
    featVect<-as.vector(featVect)

    totalMatrix[n,]<-featVect
  }

  return(totalMatrix)
}




NormalizeMBorto_CorTri<-function(maxLag=30,aaIdx,pPrim,numSeqs,lenSeqs,Trimer){
  totalMatrix<-matrix(0,nrow = numSeqs,ncol = (nrow(aaIdx)*maxLag))
  tempName<-rep(1:maxLag,nrow(aaIdx))
  tempName<-paste0("d=",tempName)
  temp2Name<-rep(row.names(aaIdx),each=maxLag)
  colnames(totalMatrix)<-paste(tempName,temp2Name)
  for(n in 1:numSeqs){
    triNuc<-Trimer[[n]]
    N<-lenSeqs[n]
    featureMartix<-matrix(ncol = nrow(aaIdx), nrow = maxLag)
    featVect=vector()
    for(d in 1:maxLag){
      vecti<-1:(N-d)
      vectj<-vecti+d
      tempMat<-(aaIdx[,triNuc[vecti]])*(aaIdx[,triNuc[vectj]])
      sigma<-apply(tempMat, 1, sum)
      numerator<-sigma*(1/(N-d))
      featVect<-rbind(featVect,numerator)
    }
    featVect<-as.vector(featVect)

    totalMatrix[n,]<-featVect
  }

  return(totalMatrix)
}


#Auto Covariance for the same property
AC_CorTri<- function(maxLag=30,aaIdx,pPrim,numSeqs,lenSeqs,Trimer){


  totalMatrix<-matrix(0,nrow = numSeqs,ncol = (nrow(aaIdx)*maxLag))
  tempName<-rep(1:maxLag,nrow(aaIdx))
  tempName<-paste0("d=",tempName)
  temp2Name<-rep(row.names(aaIdx),each=maxLag)
  colnames(totalMatrix)<-paste(tempName,temp2Name)
  for(n in 1:numSeqs){
    triNuc<-Trimer[[n]]
    N<-lenSeqs[n]
    featVect=vector()
    for(d in 1:maxLag){
      vecti<-1:(N-d)
      vectj<-vecti+d
      tempMat<-(aaIdx[,triNuc[vecti]]-pPrim[,n])*(aaIdx[,triNuc[vectj]]-pPrim[,n])
      sigma<-apply(tempMat, 1, sum)
      numerator<-sigma*(1/(N-d-1))
      featVect<-rbind(featVect,numerator)
    }
    featVect<-as.vector(featVect)

    totalMatrix[n,]<-featVect

  }
  return(totalMatrix)
}

#Cross covariance for different properties
CC_CorTri<- function(selectedProp, maxLag=30,aaIdx,pPrim,numSeqs,lenSeqs,Trimer){

  if(!is.list(selectedProp)){
    tempVect<-selectedProp
    selectedProp<-list()
    selectedProp[[1]]<-tempVect
  }
  len<-lapply(selectedProp, length)
  newProp0<-list()
  index=1
  for(i in length(len):1)
  {
    if(len[i]<2){
      selectedProp[[i]]<-NULL
      next()
    }
    if(len[i]==2){
      newVect<-c(selectedProp[[i]][2],selectedProp[[i]][1])
      newProp0[[index]]<-newVect
      index=index+1
    }
  }
  ind=1
  for (i in 1:length(selectedProp)){
    if (length(selectedProp[[ind]])>2){
      selectedProp[[ind]]<-unique(selectedProp[[ind]])
      comb<-combn(selectedProp[[ind]],2)
      newProp<-split(comb, rep(1:ncol(comb), each = nrow(comb)))
      selectedProp[ind]<-NULL

      selectedProp<-append(selectedProp,newProp)

      newProp2<-lapply(newProp, rev)

      selectedProp<-append(selectedProp,newProp2)

    }
    else {
      ind=ind+1
    }

  }
  selectedProp<-append(selectedProp,newProp0)
  selectedProp<-unique(selectedProp)

  totalMatrix<-matrix(0,nrow = numSeqs,ncol = (length(selectedProp)*maxLag))

  tempName<-rep(1:maxLag,length(selectedProp))
  tempName<-paste0("d=",tempName,"_")
  temp2Name<-rep(selectedProp,each=maxLag)
  temp2Name<-temp2Name
  colnames(totalMatrix)<-paste0(tempName,temp2Name)
  for(n in 1:numSeqs){
    triNuc<-Trimer[[n]]
    N<-lenSeqs[n]
    listComb<-list()

    featureMatrix<-vector()
    for(d in 1:maxLag){
      vecti<-1:(N-d)
      vectj<-vecti+d
      featVect<-vector(mode = "numeric",length = length(selectedProp))
      for(i in 1:length(selectedProp)){
        sigma<-sum((aaIdx[selectedProp[[i]][1],triNuc[vecti]]-pPrim[selectedProp[[i]][1],n])*
                     (aaIdx[selectedProp[[i]][2],triNuc[vectj]]-pPrim[selectedProp[[i]][2],n]))
        featVect[i]<-sigma
      }
      featVect<-featVect/(N-d-1)
      featureMatrix<-rbind(featureMatrix,featVect)

    }
    featureMatrix<-as.vector(featureMatrix)
    totalMatrix[n,]<-featureMatrix

  }
  return(totalMatrix)
}


ACC_CorTri<- function(selectedProp, maxLag=30,aaIdx,pPrim,numSeqs,lenSeqs,Trimer){


  if(!is.list(selectedProp)){
    tempVect<-selectedProp
    selectedProp<-list()
    selectedProp[[1]]<-tempVect
  }
  len<-lapply(selectedProp, length)
  newProp0<-list()
  index=1
  lens<-length(len)
  for(i in lens:1)
  {

    if(len[i]==1){
      newVect<-c(selectedProp[[i]],selectedProp[[i]])
      newProp0[[index]]<-newVect
      index=index+1
      selectedProp[i]<-NULL
      next()
    }

    if(len[i]==0){
      selectedProp[i]<-NULL
      next()
    }

    if(len[i]==2){
      newVect<-c(selectedProp[[i]][2],selectedProp[[i]][1])
      newProp0[[index]]<-newVect
      index=index+1
    }
  }
  ind=1
  for (i in 1:length(selectedProp)){
    if (length(selectedProp[[ind]])>2){
      vect<-selectedProp[[ind]]
      names(vect)<-names(selectedProp[[ind]])
      vect<-unique(vect)
      comb<-combn(vect,2)
      newProp<-split(comb, rep(1:ncol(comb), each = nrow(comb)))
      selectedProp[ind]<-NULL
      selectedProp<-append(selectedProp,newProp)

      newProp2<-lapply(newProp, rev)

      selectedProp<-append(selectedProp,newProp2)

    } else {
      ind=ind+1
    }
  }
  selectedProp<-unique(selectedProp)
  unlistProp<-unique(unlist(selectedProp))
  for(i in 1:length(unlistProp)){
    newProp<-c(unlistProp[i],unlistProp[i])
    n<-length(selectedProp)
    selectedProp[[(n+1)]]<-newProp
  }
  selectedProp<-append(selectedProp,newProp0)
  selectedProp<-unique(selectedProp)

  totalMatrix<-matrix(0,nrow = numSeqs,ncol = (length(selectedProp)*maxLag))
  tempName<-rep(1:maxLag,length(selectedProp))
  tempName<-paste0("d=",tempName,"_")
  temp2Name<-rep(selectedProp,each=maxLag)
  temp2Name<-temp2Name
  colnames(totalMatrix)<-paste0(tempName,temp2Name)
  for(n in 1:numSeqs){
    triNuc<-Trimer[[n]]
    N<-lenSeqs[n]

    featureMatrix<-vector()
    for(d in 1:maxLag){
      vecti<-1:(N-d)
      vectj<-vecti+d
      featVect<-vector(mode = "numeric",length = length(selectedProp))
      for(i in 1:length(selectedProp)){
        sigma<-sum((aaIdx[selectedProp[[i]][1],triNuc[vecti]]-pPrim[selectedProp[[i]][1],n])*
                     (aaIdx[selectedProp[[i]][2],triNuc[vectj]]-pPrim[selectedProp[[i]][2],n]))
        featVect[i]<-sigma

      }
      featVect<-featVect/(N-d-1)
      featureMatrix<-rbind(featureMatrix,featVect)

    }
    featureMatrix<-as.vector(featureMatrix)
    totalMatrix[n,]<-featureMatrix

  }

  return(totalMatrix)

}
