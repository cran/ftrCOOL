#' Amino Acid Autocorrelation-Autocovariance (AAutoCor)
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
#' @param seqs is a FASTA file with amino acid sequences. Each sequence starts
#' with a '>' character. Also, seqs could be a string vector. Each element of the vector is a peptide/protein sequence.
#'
#'
#' @param selectedAAidx Function takes as input the physicochemical properties. Users select the properties by their ids
#' or indices in the aaIndex2 file. This parameter could be a vector or a list of amino acid indices.
#' The default values of the vector are the 'CIDH920105','BHAR880101','CHAM820101','CHAM820102','CHOC760101','BIGC670101','CHAM810101','DAYM780201'
#' ids in the aaIndex2 file.
#'
#'
#' @param maxlag This parameter shows the maximum gap between two amino acids. The gaps change from 1 to maxlag (the maximum lag).
#'
#' @param type could be 'Moran', 'Greay', 'NormalizeMBorto', 'AC', 'CC', or 'ACC'. Also, it could be any combination of them.
#'
#' @param threshold is a number between (0 , 1]. In selectedAAidx, indices with a correlation higher than the threshold will be deleted.
#' The default value is 1.
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @return This function returns a feature matrix. The number of columns in the matrix changes depending on the chosen autocorrelation or autocovariance types and nlag parameter.
#' The output is a matrix. The number of rows shows the number of sequences.
#'
#' @export
#'
#' @examples
#'
#' filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
#' mat1<-AAutoCor(seqs=filePrs,maxlag=20,threshold=0.9,
#' type=c("Moran","Geary","NormalizeMBorto","AC"))
#'
#' mat2<-AAutoCor(seqs=filePrs,maxlag=20,threshold=0.9,selectedAAidx=
#' list(c('CIDH920105','BHAR880101','CHAM820101','CHAM820102'),c('CHOC760101','BIGC670101')
#' ,c('CHAM810101','DAYM780201')),type=c("AC","CC","ACC"))
#'


AAutoCor<-function(seqs,selectedAAidx=list(c('CIDH920105','BHAR880101','CHAM820101','CHAM820102','CHOC760101','BIGC670101','CHAM810101','DAYM780201')),maxlag=3,threshold=1,type=c("Moran","Geary","NormalizeMBorto","AC","CC","ACC"),label=c()){


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
  lenSeqs<-sapply(seqs, nchar)
  charseq<-lapply(seqs, function(i) unlist(strsplit(i,split = "")))

  aaIdxAd<-paste0(path.pack,"/standardizedAAindex_555.csv")
  aaIdx<-read.csv(aaIdxAd)


  rNumAAIdx<-unlist(selectedAAidx)
  unqrNumAAIdx<-unique(rNumAAIdx)
  row.names(aaIdx)<-aaIdx[,1]
  aaIdx<-as.matrix(aaIdx)


  selectedAAidxNames=list()
  if(is.list(selectedAAidx))
  {
    for(i in 1:length(selectedAAidx)){
      if(is.numeric(selectedAAidx[[i]])){
        names(selectedAAidx[[i]])<-row.names(aaIdx)[selectedAAidx[[i]]]
      }
      else if(is.character(selectedAAidx[[i]])){
        names(selectedAAidx[[i]])<-selectedAAidx[[i]]
      } else{
        stop('ERROR: selectedAAidx should be Ids or indices of the selected properties in the aaIndex2 file')
      }
    }

    selectedAAidxNames<-lapply(selectedAAidx, function(x){
      names(x)
    })
  }
  else if(is.vector(selectedAAidx)){
    if(is.numeric(selectedAAidx)){
      names(selectedAAidx)<-row.names(aaIdx)[selectedAAidx]
    } else if(is.character(selectedAAidx)){
      names(selectedAAidx)<-selectedAAidx
    } else{
      stop('ERROR: selectedAAidx should be Ids or indices of the selected properties in the aaIndex2 file')
    }
    selectedAAidxNames=c()
    selectedAAidxNames<-names(selectedAAidx)
  }
  else
    stop("ERROR: selectedAAidx should be a vector or a list")


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
      selectedAAidxNames<-lapply(selectedAAidxNames, function(x){
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
      Pprim[j,i]<-sum(aaIdx[i,charseq[[j]]])/lenSeqs[j]
    }
  }

  Pprim<-t(Pprim)
  row.names(Pprim)<-row.names(aaIdx)
  totalfeatureMatrix=vector()
  if("Moran" %in% type){
    featureMatrix<-Moran_corAA(maxLag = maxlag,aaIdx = aaIdx,pPrim = Pprim,numSeqs = numSeqs, lenSeqs=lenSeqs,charact = charseq)
    colnames(featureMatrix)<-paste("Mo",colnames(featureMatrix))
    totalfeatureMatrix<-cbind(totalfeatureMatrix,featureMatrix)
  }
  if("Geary" %in% type){
    featureMatrix<-Geary_corAA(maxLag = maxlag,aaIdx = aaIdx,pPrim = Pprim,numSeqs = numSeqs,lenSeqs = lenSeqs,charact = charseq)
    colnames(featureMatrix)<-paste("Ge",colnames(featureMatrix))
    totalfeatureMatrix<-cbind(totalfeatureMatrix,featureMatrix)
  }
  if("NormalizeMBorto" %in% type){
    featureMatrix<-NormalizeMBorto_corAA(maxLag = maxlag,aaIdx = aaIdx,pPrim = Pprim,numSeqs = numSeqs,lenSeqs = lenSeqs,charact = charseq)
    colnames(featureMatrix)<-paste("No",colnames(featureMatrix))
    totalfeatureMatrix<-cbind(totalfeatureMatrix,featureMatrix)
  }
  if("AC" %in% type){
    featureMatrix<-AC_CorAA(maxLag = maxlag,aaIdx = aaIdx,pPrim = Pprim,numSeqs = numSeqs,lenSeqs = lenSeqs,charact = charseq)
    colnames(featureMatrix)<-paste("AC",colnames(featureMatrix))
    totalfeatureMatrix<-cbind(totalfeatureMatrix,featureMatrix)
  }
  if("CC" %in% type){
    featureMatrix<-CC_CorAA(selectedProp = selectedAAidxNames ,maxLag = maxlag,aaIdx = aaIdx,pPrim = Pprim,numSeqs = numSeqs,lenSeqs = lenSeqs,charact = charseq)
    colnames(featureMatrix)<-paste("CC",colnames(featureMatrix))
    totalfeatureMatrix<-cbind(totalfeatureMatrix,featureMatrix)
  }
  if("ACC" %in% type){
    featureMatrix<-ACC_CorAA(selectedProp = selectedAAidxNames ,maxLag = maxlag,aaIdx = aaIdx,pPrim = Pprim,numSeqs = numSeqs,lenSeqs = lenSeqs,charact = charseq)
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

Moran_corAA<-function(maxLag=30,aaIdx,pPrim,numSeqs,lenSeqs,charact){

  totalMatrix<-matrix(0,nrow = numSeqs,ncol = (nrow(aaIdx)*maxLag))
  tempName<-rep(1:maxLag,nrow(aaIdx))
  tempName<-paste0("d=",tempName)
  temp2Name<-rep(row.names(aaIdx),each=maxLag)
  colnames(totalMatrix)<-paste(tempName,temp2Name)
  for(n in 1:numSeqs){
     aa<-charact[[n]]
     N<-lenSeqs[n]
     featureMartix<-matrix(ncol = nrow(aaIdx), nrow = maxLag)
     tempMat<-(aaIdx[,aa]-pPrim[,n])^2
     sigmaRow<-apply(tempMat,1,sum)
     denominator<-(1/lenSeqs[n])*sigmaRow
     featVect=vector()
     for(d in 1:maxLag){
        vecti<-1:(N-d)
        vectj<-vecti+d
        tempMat<-(aaIdx[,aa[vecti]]-pPrim[,n])*(aaIdx[,aa[vectj]]-pPrim[,n])
        sigma<-apply(tempMat, 1, sum)
        numerator<-sigma*(1/(N-d))/denominator
        featVect<-rbind(featVect,numerator)
     }
     featVect<-as.vector(featVect)

     totalMatrix[n,]<-featVect

  }
  return(totalMatrix)
}

Geary_corAA<-function(maxLag=30,aaIdx,pPrim,numSeqs,lenSeqs,charact){
  totalMatrix<-matrix(0,nrow = numSeqs,ncol = (nrow(aaIdx)*maxLag))
  tempName<-rep(1:maxLag,nrow(aaIdx))
  tempName<-paste0("d=",tempName)
  temp2Name<-rep(row.names(aaIdx),each=maxLag)
  colnames(totalMatrix)<-paste(tempName,temp2Name)
  for(n in 1:numSeqs){
    aa<-charact[[n]]
    N<-lenSeqs[n]
    featureMartix<-matrix(ncol = nrow(aaIdx), nrow = maxLag)
    tempMat<-(aaIdx[,aa]-pPrim[,n])^2
    sigmaRow<-apply(tempMat,1,sum)
    denominator<-(1/lenSeqs[n])*sigmaRow
    featVect=vector()
    for(d in 1:maxLag){
      vecti<-1:(N-d)
      vectj<-vecti+d
      tempMat<-(aaIdx[,aa[vecti]]-aaIdx[,aa[vectj]])^2
      sigma<-apply(tempMat, 1, sum)
      numerator<-sigma*(1/(N-d))/denominator
      featVect<-rbind(featVect,numerator)
    }
    featVect<-as.vector(featVect)

    totalMatrix[n,]<-featVect
  }

  return(totalMatrix)
}




NormalizeMBorto_corAA<-function(maxLag=30,aaIdx,pPrim,numSeqs,lenSeqs,charact){
  totalMatrix<-matrix(0,nrow = numSeqs,ncol = (nrow(aaIdx)*maxLag))
  tempName<-rep(1:maxLag,nrow(aaIdx))
  tempName<-paste0("d=",tempName)
  temp2Name<-rep(row.names(aaIdx),each=maxLag)
  colnames(totalMatrix)<-paste(tempName,temp2Name)
  for(n in 1:numSeqs){
    aa<-charact[[n]]
    N<-lenSeqs[n]
    featureMartix<-matrix(ncol = nrow(aaIdx), nrow = maxLag)
    featVect=vector()
    for(d in 1:maxLag){
      vecti<-1:(N-d)
      vectj<-vecti+d
      tempMat<-(aaIdx[,aa[vecti]])*(aaIdx[,aa[vectj]])
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
AC_CorAA<- function(maxLag=30,aaIdx,pPrim,numSeqs,lenSeqs,charact){


  totalMatrix<-matrix(0,nrow = numSeqs,ncol = (nrow(aaIdx)*maxLag))
  tempName<-rep(1:maxLag,nrow(aaIdx))
  tempName<-paste0("d=",tempName)
  temp2Name<-rep(row.names(aaIdx),each=maxLag)
  colnames(totalMatrix)<-paste(tempName,temp2Name)
  for(n in 1:numSeqs){
    aa<-charact[[n]]
    N<-lenSeqs[n]
    featVect=vector()
    for(d in 1:maxLag){
      vecti<-1:(N-d)
      vectj<-vecti+d
      tempMat<-(aaIdx[,aa[vecti]]-pPrim[,n])*(aaIdx[,aa[vectj]]-pPrim[,n])
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
CC_CorAA<- function(selectedProp, maxLag=30,aaIdx,pPrim,numSeqs,lenSeqs,charact){

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
    aa<-charact[[n]]
    N<-lenSeqs[n]
    listComb<-list()

    featureMatrix<-vector()
    for(d in 1:maxLag){
      vecti<-1:(N-d)
      vectj<-vecti+d
      featVect<-vector(mode = "numeric",length = length(selectedProp))
      for(i in 1:length(selectedProp)){
        sigma<-sum((aaIdx[selectedProp[[i]][1],aa[vecti]]-pPrim[selectedProp[[i]][1],n])*
          (aaIdx[selectedProp[[i]][2],aa[vectj]]-pPrim[selectedProp[[i]][2],n]))
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


ACC_CorAA<- function(selectedProp, maxLag=30,aaIdx,pPrim,numSeqs,lenSeqs,charact){


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
    aa<-charact[[n]]
    N<-lenSeqs[n]

    featureMatrix<-vector()
    for(d in 1:maxLag){
      vecti<-1:(N-d)
      vectj<-vecti+d
      featVect<-vector(mode = "numeric",length = length(selectedProp))
      for(i in 1:length(selectedProp)){
        sigma<-sum((aaIdx[selectedProp[[i]][1],aa[vecti]]-pPrim[selectedProp[[i]][1],n])*
                     (aaIdx[selectedProp[[i]][2],aa[vectj]]-pPrim[selectedProp[[i]][2],n]))
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
