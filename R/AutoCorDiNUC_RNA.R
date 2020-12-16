#' Di riboNucleotide Autocorrelation-Autocovariance (AutoCorDiNUC_RNA)
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
#' @param seqs is a FASTA file containing ribonucleic acid(RNA) sequences. The sequences start
#' with '>'. Also, seqs could be a string vector. Each element of the vector is a RNA sequence.
#'
#'
#' @param selectedIdx function takes as input the physicochemical properties. Users select the properties by their ids
#' or indices in the DI_RNA file. This parameter could be a vector or a list of di-ribonucleic acid indices.
#' The default value of this parameter is a vector with ("Rise (RNA)", "Roll (RNA)", "Shift (RNA)", "Slide (RNA)", "Tilt (RNA)","Twist (RNA)") ids.
#'
#'
#' @param maxlag This parameter shows the maximum gap between two amino acids. The gaps change from 1 to maxlag (the maximum lag).
#'
#' @param type could be 'Moran', 'Greay', 'NormalizeMBorto', 'AC', 'CC', or 'ACC'. Also, it could be any combination of them.
#'
#' @param threshold is a number between (0 to 1]. In selectedIdx, indices with a correlation higher than the threshold will be deleted.The default value is 1.
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#'
#' @return This function returns a feature matrix. The number of columns in the matrix changes depending on the chosen autocorrelation or autocovariance types and nlag parameter.
#' The output is a matrix. The number of rows shows the number of sequences.
#'
#' @export
#'
#' @examples
#'
#' fileLNC<-system.file("extdata/Carica_papaya101RNA.txt",package="ftrCOOL")
#' fileLNC<-fa.read(fileLNC,alphabet="rna")
#' fileLNC<-fileLNC[1:20]
#' mat1<-AutoCorDiNUC_RNA(seqs=fileLNC,maxlag=20,type=c("Moran","Geary"))
#'

AutoCorDiNUC_RNA<-function(seqs,selectedIdx=c("Rise (RNA)", "Roll (RNA)", "Shift (RNA)", "Slide (RNA)", "Tilt (RNA)","Twist (RNA)"),maxlag=3,threshold=1,type=c("Moran","Geary","NormalizeMBorto","AC","CC","ACC"),label=c()){


  path.pack=system.file("extdata",package="ftrCOOL")
  if(length(seqs)==1&&file.exists(seqs)){
    seqs<-fa.read(seqs,alphabet="rna")
    seqs_Lab<-alphabetCheck(seqs,alphabet = "rna",label)

    seqs<-seqs_Lab[[1]]
    label<-seqs_Lab[[2]]
  }
  else if(is.vector(seqs)){
    seqs<-sapply(seqs,toupper)

    seqs_Lab<-alphabetCheck(seqs,alphabet = "rna",label)

    seqs<-seqs_Lab[[1]]
    label<-seqs_Lab[[2]]

  }
  else {
    stop("ERROR: Input sequence is not in the correct format. It should be a FASTA file or a string vector.")
  }


  numSeqs<-length(seqs)
  lenSeqs<-sapply(seqs, nchar)
  lenSeqs<-lenSeqs-1
  Di<-lapply(seqs, function(i) {
    chars<-unlist(strsplit(i,split = ""))
    N=nchar(i)
    temp1<-chars[1:(N-1)]
    temp2<-chars[2:N]
    dimers<-paste0(temp1,temp2)
    return(dimers)
  })


  #aaIdxAd<-paste0(path.pack,"/standardizeNUCidx2.csv")
  #aaIdxAd<-paste0(path.pack,"/DI_DNA.csv")
  aaIdxAd<-paste0(path.pack,"/DI_RNA.csv")
  aaIdx<-read.csv(aaIdxAd)


  rNumAAIdx<-unlist(selectedIdx)
  unqrNumAAIdx<-unique(rNumAAIdx)
  row.names(aaIdx)<-aaIdx[,1]
  aaIdx<-as.matrix(aaIdx)

  selectedIdxNames=list()
  if(is.list(selectedIdx))
  {
    for(i in 1:length(selectedIdx)){
      if(is.numeric(selectedIdx[[i]])){
        names(selectedIdx[[i]])<-row.names(aaIdx)[selectedIdx[[i]]]
      }
      else if(is.character(selectedIdx[[i]])){
        names(selectedIdx[[i]])<-selectedIdx[[i]]
      } else{
        stop('ERROR: selectedIdx should be Ids or indices of the selected properties in the di-ribonucleotide index file')
      }
    }

    selectedIdxNames<-lapply(selectedIdx, function(x){
      names(x)
    })
  }
  else if(is.vector(selectedIdx)){
    if(is.numeric(selectedIdx)){
      names(selectedIdx)<-row.names(aaIdx)[selectedIdx]
    } else if(is.character(selectedIdx)){
      names(selectedIdx)<-selectedIdx
    } else{
      stop('ERROR: selectedIdx should be Ids or indices of the selected properties in the di-ribonucleotide index file')
    }
    selectedIdxNames=c()
    selectedIdxNames<-names(selectedIdx)
  }
  else{
    stop("ERROR: selectedIdx should be a vector or a list")
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
      selectedIdxNames<-lapply(selectedIdxNames, function(x){
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



  ## average all amino acids in seq for each property
  ##is a vector with length property number
  Pprim <- matrix(0, ncol = numAAidx,nrow = numSeqs)


  # Is a matrix each row belong to a sequence each col belong to each property
  for( j in 1:numSeqs){
    for(i in 1:numAAidx)
    {
      Pprim[j,i]<-sum(aaIdx[i,Di[[j]]])/lenSeqs[j]
    }
  }

  Pprim<-t(Pprim)
  row.names(Pprim)<-row.names(aaIdx)
  totalfeatureMatrix=vector()
  if("Moran" %in% type){
    featureMatrix<-Moran_CorDi_RNA(maxLag = maxlag,aaIdx = aaIdx,pPrim = Pprim,numSeqs=numSeqs,lenSeqs=lenSeqs,Dimer = Di)
    colnames(featureMatrix)<-paste("Mo",colnames(featureMatrix))
    totalfeatureMatrix<-cbind(totalfeatureMatrix,featureMatrix)
  }
  if("Geary" %in% type){
    featureMatrix<-Geary_CorDi_RNA(maxLag = maxlag,aaIdx = aaIdx,pPrim = Pprim,numSeqs=numSeqs,lenSeqs=lenSeqs,Dimer = Di)
    colnames(featureMatrix)<-paste("Ge",colnames(featureMatrix))
    totalfeatureMatrix<-cbind(totalfeatureMatrix,featureMatrix)
  }
  if("NormalizeMBorto" %in% type){
    featureMatrix<-NormalizeMBorto_CorDi_RNA(maxLag = maxlag,aaIdx = aaIdx,pPrim = Pprim,numSeqs=numSeqs,lenSeqs=lenSeqs,Dimer = Di)
    colnames(featureMatrix)<-paste("No",colnames(featureMatrix))
    totalfeatureMatrix<-cbind(totalfeatureMatrix,featureMatrix)
  }
  if("AC" %in% type){
    featureMatrix<-AC_CorDi_RNA(maxLag = maxlag,aaIdx = aaIdx,pPrim = Pprim,numSeqs=numSeqs,lenSeqs=lenSeqs,Dimer = Di)
    colnames(featureMatrix)<-paste("AC",colnames(featureMatrix))
    totalfeatureMatrix<-cbind(totalfeatureMatrix,featureMatrix)
  }
  if("CC" %in% type){
    featureMatrix<-CC_CorDi_RNA(selectedProp = selectedIdxNames ,maxLag = maxlag,aaIdx = aaIdx,pPrim = Pprim,numSeqs=numSeqs,lenSeqs=lenSeqs,Dimer = Di)
    colnames(featureMatrix)<-paste("CC",colnames(featureMatrix))
    totalfeatureMatrix<-cbind(totalfeatureMatrix,featureMatrix)
  }
  if("ACC" %in% type){
    featureMatrix<-ACC_CorDi_RNA(selectedProp = selectedIdxNames ,maxLag = maxlag,aaIdx = aaIdx,pPrim = Pprim,numSeqs=numSeqs,lenSeqs=lenSeqs,Dimer = Di)
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

Moran_CorDi_RNA<-function(maxLag=30,aaIdx,pPrim,numSeqs,lenSeqs,Dimer){

  totalMatrix<-matrix(0,nrow = numSeqs,ncol = (nrow(aaIdx)*maxLag))
  tempName<-rep(1:maxLag,nrow(aaIdx))
  tempName<-paste0("d=",tempName)
  temp2Name<-rep(row.names(aaIdx),each=maxLag)
  colnames(totalMatrix)<-paste(tempName,temp2Name)
  for(n in 1:numSeqs){
    diNuc<-Dimer[[n]]
    N<-lenSeqs[n]
    featureMartix<-matrix(ncol = nrow(aaIdx), nrow = maxLag)
    tempMat<-(aaIdx[,diNuc]-pPrim[,n])^2
    sigmaRow<-apply(tempMat,1,sum)
    denominator<-(1/lenSeqs[n])*sigmaRow
    featVect=vector()
    for(d in 1:maxLag){
      vecti<-1:(N-d)
      vectj<-vecti+d
      tempMat<-(aaIdx[,diNuc[vecti]]-pPrim[,n])*(aaIdx[,diNuc[vectj]]-pPrim[,n])
      sigma<-apply(tempMat, 1, sum)
      numerator<-sigma*(1/(N-d))/denominator
      featVect<-rbind(featVect,numerator)
    }
    featVect<-as.vector(featVect)

    totalMatrix[n,]<-featVect

  }
  return(totalMatrix)
}

Geary_CorDi_RNA<-function(maxLag=30,aaIdx,pPrim,numSeqs,lenSeqs,Dimer){
  totalMatrix<-matrix(0,nrow = numSeqs,ncol = (nrow(aaIdx)*maxLag))
  tempName<-rep(1:maxLag,nrow(aaIdx))
  tempName<-paste0("d=",tempName)
  temp2Name<-rep(row.names(aaIdx),each=maxLag)
  colnames(totalMatrix)<-paste(tempName,temp2Name)
  for(n in 1:numSeqs){
    diNuc<-Dimer[[n]]
    N<-lenSeqs[n]
    featureMartix<-matrix(ncol = nrow(aaIdx), nrow = maxLag)
    tempMat<-(aaIdx[,diNuc]-pPrim[,n])^2
    sigmaRow<-apply(tempMat,1,sum)
    denominator<-(1/lenSeqs[n])*sigmaRow
    featVect=vector()
    for(d in 1:maxLag){
      vecti<-1:(N-d)
      vectj<-vecti+d
      tempMat<-(aaIdx[,diNuc[vecti]]-aaIdx[,diNuc[vectj]])^2
      sigma<-apply(tempMat, 1, sum)
      numerator<-sigma*(1/(N-d))/denominator
      featVect<-rbind(featVect,numerator)
    }
    featVect<-as.vector(featVect)

    totalMatrix[n,]<-featVect
  }

  return(totalMatrix)
}




NormalizeMBorto_CorDi_RNA<-function(maxLag=30,aaIdx,pPrim,numSeqs,lenSeqs,Dimer){
  totalMatrix<-matrix(0,nrow = numSeqs,ncol = (nrow(aaIdx)*maxLag))
  tempName<-rep(1:maxLag,nrow(aaIdx))
  tempName<-paste0("d=",tempName)
  temp2Name<-rep(row.names(aaIdx),each=maxLag)
  colnames(totalMatrix)<-paste(tempName,temp2Name)
  for(n in 1:numSeqs){
    diNuc<-Dimer[[n]]
    N<-lenSeqs[n]
    featureMartix<-matrix(ncol = nrow(aaIdx), nrow = maxLag)
    featVect=vector()
    for(d in 1:maxLag){
      vecti<-1:(N-d)
      vectj<-vecti+d
      tempMat<-(aaIdx[,diNuc[vecti]])*(aaIdx[,diNuc[vectj]])
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
AC_CorDi_RNA<- function(maxLag=30,aaIdx,pPrim,numSeqs,lenSeqs,Dimer){


  totalMatrix<-matrix(0,nrow = numSeqs,ncol = (nrow(aaIdx)*maxLag))
  tempName<-rep(1:maxLag,nrow(aaIdx))
  tempName<-paste0("d=",tempName)
  temp2Name<-rep(row.names(aaIdx),each=maxLag)
  colnames(totalMatrix)<-paste(tempName,temp2Name)
  for(n in 1:numSeqs){
    diNuc<-Dimer[[n]]
    N<-lenSeqs[n]
    featVect=vector()
    for(d in 1:maxLag){
      vecti<-1:(N-d)
      vectj<-vecti+d
      tempMat<-(aaIdx[,diNuc[vecti]]-pPrim[,n])*(aaIdx[,diNuc[vectj]]-pPrim[,n])
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
CC_CorDi_RNA<- function(selectedProp, maxLag=30,aaIdx,pPrim,numSeqs,lenSeqs,Dimer){

  if(!is.list(selectedProp)){
    tempVect<-selectedProp
    selectedProp<-list()
    selectedProp[[1]]<-tempVect
  }
  len<-lapply(selectedProp, length)
  for(i in length(len):1)
  {
    if(len[i]<2){
      selectedProp[[i]]<-NULL
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
  selectedProp<-unique(selectedProp)

  totalMatrix<-matrix(0,nrow = numSeqs,ncol = (length(selectedProp)*maxLag))

  tempName<-rep(1:maxLag,length(selectedProp))
  tempName<-paste0("d=",tempName,"_")
  temp2Name<-rep(selectedProp,each=maxLag)
  temp2Name<-temp2Name
  colnames(totalMatrix)<-paste0(tempName,temp2Name)
  for(n in 1:numSeqs){
    diNuc<-Dimer[[n]]
    N<-lenSeqs[n]
    listComb<-list()

    featureMatrix<-vector()
    for(d in 1:maxLag){
      vecti<-1:(N-d)
      vectj<-vecti+d
      featVect<-vector(mode = "numeric",length = length(selectedProp))
      for(i in 1:length(selectedProp)){
        sigma<-sum((aaIdx[selectedProp[[i]][1],diNuc[vecti]]-pPrim[selectedProp[[i]][1],n])*
                     (aaIdx[selectedProp[[i]][2],diNuc[vectj]]-pPrim[selectedProp[[i]][2],n]))
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


ACC_CorDi_RNA<- function(selectedProp, maxLag=30,aaIdx,pPrim,numSeqs,lenSeqs,Dimer){


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
    #print(len[i])
    if(len[i]==0){
      #print(len[i])
      selectedProp[i]<-NULL
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
    diNuc<-Dimer[[n]]
    N<-lenSeqs[n]

    featureMatrix<-vector()
    for(d in 1:maxLag){
      vecti<-1:(N-d)
      vectj<-vecti+d
      featVect<-vector(mode = "numeric",length = length(selectedProp))
      for(i in 1:length(selectedProp)){
        sigma<-sum((aaIdx[selectedProp[[i]][1],diNuc[vecti]]-pPrim[selectedProp[[i]][1],n])*
                     (aaIdx[selectedProp[[i]][2],diNuc[vectj]]-pPrim[selectedProp[[i]][2],n]))
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
