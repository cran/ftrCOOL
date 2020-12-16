#' Composition_Transition_Distribution (CTD)
#'
#' This function calculates the composition, transition, and distribution for each sequence.
#'
#' @references Dubchak, Inna, et al. "Prediction of protein folding class using global description of amino acid sequence." Proceedings of the National Academy of Sciences 92.19 (1995): 8700-8704.
#'
#' @param seqs is a FASTA file with amino acid sequences. Each sequence starts
#' with a '>' character. Also, seqs could be a string vector. Each element of the vector is a peptide/protein sequence.
#'
#' @param normalized is a logical parameter. When it is FALSE, the return value of the function does not change. Otherwise, the return value is normalized using the length of the sequence.
#'
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#'
#' @return Output is a combination of three different matrices: Composition, Transition, and Distribution.
#' You can obtain any of the three matrices by executing the corresponding function, i.e., \link{CTDC}, \link{CTDT}, and \link{CTDD}.
#'
#'
#'
#' @export
#'
#' @examples
#'
#' filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
#' CTDtotal<-CTD(seqs=filePrs,normalized=FALSE)
#'

CTD<-function(seqs,normalized=FALSE,label=c()){




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

  }else {
    stop("ERROR: Input sequence is not in the correct format. It should be a FASTA file or a string vector.")
  }




  numSeqs<-length(seqs)
  numGrp=7

  featureMatrix<-matrix(0,ncol = ((3*numGrp)+(3*numGrp)+(15*numGrp)),nrow = numSeqs)

  temp1<-paste0("P",rep(1:numGrp,each=3))
  temp2<-paste0(".G",rep(1:3,numGrp))
  namesC<-paste0("C(",temp1,temp2,")")


  temp1NaT<-paste0("P",rep(1:numGrp,each=3))
  temp2NaT<-rep(c("(1,2)","(1,3)","(2,3)"),7)
  namesT<-paste0("T(",temp1NaT,temp2NaT,")")

  nameP1<-paste("P",rep(1:numGrp,each=(3*5)),sep = "")
  nameP2<-paste(".G",rep(c(1,2,3),each=5,numGrp),sep="")
  nameP3<-paste(rep(c(".First",".25%",".50%",".75%",".100%"),(numGrp*3)))
  tempnName<-paste0(nameP1,nameP2,nameP3)
  namesD<-paste0("D(",tempnName,")")

  names<-c(namesC,namesT,namesD)
  colnames(featureMatrix)<-names


  group1<-list("hydrophobicity_PRAM900101"= c("R","K","E","D","Q","N"),
               "normwaalsvolume"= c("G","A","S","T","P","D","C"),
               "polarity"=        c("L","I","F","W","C","M","V","Y"),
               "polarizability"=  c("G","A","S","D","T"),
               "charge"=          c("K","R"),
               "secondarystruct"= c("E","A","L","M","Q","K","R","H"),
               "solventaccess"=   c("A","L","F","C","G","I","V","W"))

  group2<-list("hydrophobicity_PRAM900101"= c("G","A","S","T","P","H","Y"),
               "normwaalsvolume"= c("N","V","E","Q","I","L"),
               "polarity"=        c("P","A","T","G","S"),
               "polarizability"=  c("C","P","N","V","E","Q","I","L"),
               "charge"=          c("A","N","C","Q","G","H","I","L","M","F","P","S","T","W","Y","V"),
               "secondarystruct"= c("V","I","Y","C","W","F","T"),
               "solventaccess"=  c("R","K","Q","E","N","D"))

  group3<-list("hydrophobicity_PRAM900101"= c("C","L","V","I","M","F","W"),
               "normwaalsvolume"= c("M","H","K","F","R","Y","W"),
               "polarity"=        c("H","Q","R","K","N","E","D"),
               "polarizability"=  c("K","M","H","F","R","Y","W"),
               "charge"=          c("D","E"),
               "secondarystruct"= c("G","N","P","S","D"),
               "solventaccess"=   c("M","S","P","T","H","Y"))
  groups<-list("grp1"=group1,"grp2"=group2,"grp3"=group3)

  properties<-c("hydrophobicity_PRAM900101", "normwaalsvolume",
                "polarity", "polarizability", "charge", "secondarystruct", "solventaccess")


  numGrp=7


  stC<-1
  enC<-numGrp*3

  stT<-enC+1
  enT<-numGrp*6

  stD<-enT+1
  enD<-ncol(featureMatrix)

  for(n in 1:numSeqs){

    seq=seqs[n]

    aa=unlist(strsplit(seq,split = ""))
    lenaa<-length(aa)

    g1 <- lapply(group1, function(g) which(aa %in% g))

    g2 <- lapply(group2, function(g) which(aa %in% g))

    g3 <- lapply(group3, function(g) which(aa %in% g))

    CTDCvectg1<-lapply(g1, length)
    CTDCvectg2<-lapply(g2, length)
    CTDCvectg3<-lapply(g3, length)
    ctdC<-mapply(c,CTDCvectg1,CTDCvectg2,CTDCvectg3, SIMPLIFY=FALSE)
    ctdC<-unlist(ctdC)


    G <- vector("list", numGrp)

    for (i in 1:numGrp) G[[i]] <- rep(NA, lenaa)
    for (i in 1:numGrp) {
      try(G[[i]][which(aa %in% group1[[i]])] <- "1")
      try(G[[i]][which(aa %in% group2[[i]])] <- "2")
      try(G[[i]][which(aa %in% group3[[i]])] <- "3")
    }
    for (i in 1:numGrp) G[[i]] <- paste(G[[i]][-lenaa], G[[i]][-1], sep = "")
    mat<-matrix(0,ncol = 9,nrow = numGrp)
    colnames(mat)<-c("11","12","13","21","22","23","31","32","33")

    tabG<-lapply(G,table)

    for(i in 1:numGrp){
      mat[i,names(tabG[[i]])]=tabG[[i]]
    }

    mtx<-matrix(0,nrow = numGrp,ncol = 3)
    colnames(mtx)<-c("1,2","1,3","2,3")
    for(i in 1:numGrp){
      mtx[i,"1,2"]<-mat[i,2]+mat[i,4]
      mtx[i,"1,3"]<-mat[i,3]+mat[i,7]
      mtx[i,"2,3"]<-mat[i,6]+mat[i,8]
    }
    mtx<-t(mtx)
    ctdT<-as.vector(mtx)
    if(normalized==TRUE){
      ctdT<-ctdT/lenaa
    }



    featureMatrix[n,stT:enT]<-ctdT

    index.0<-which(ctdC==0)



    ctdF=rep(1,length(ctdC))
    ctdC1.4=round((1/4)*ctdC)
    ctdC2.4=round((1/2)*ctdC)
    ctdC3.4=round((3/4)*ctdC)

    ctdC[index.0]<-NA
    ctdF[index.0]<-NA
    ctdC1.4[ctdC1.4==0]<-NA
    ctdC2.4[ctdC2.4==0]<-NA
    ctdC3.4[ctdC3.4==0]<-NA




    ctdD<-c()

    for(i in 1:numGrp){
      v<-vector(mode = "numeric",length = 15)
      ind<-((i-1)*3)+1

      v[1]=g1[[i]][ctdF[ind]]
      v[2]=g1[[i]][ctdC1.4[ind]]
      v[3]=g1[[i]][ctdC2.4[ind]]
      v[4]=g1[[i]][ctdC3.4[ind]]
      v[5]=g1[[i]][ctdC[ind]]

      v[6]=g2[[i]][ctdF[(ind+1)]]
      v[7]=g2[[i]][ctdC1.4[(ind+1)]]
      v[8]=g2[[i]][ctdC2.4[(ind+1)]]
      v[9]=g2[[i]][ctdC3.4[(ind+1)]]
      v[10]=g2[[i]][ctdC[(ind+1)]]


      v[11]=g3[[i]][ctdF[(ind+2)]]
      v[12]=g3[[i]][ctdC1.4[(ind+2)]]
      v[13]=g3[[i]][ctdC2.4[(ind+2)]]
      v[14]=g3[[i]][ctdC3.4[(ind+2)]]
      v[15]=g3[[i]][ctdC[(ind+2)]]


      v[is.na(v)]<-0
      ctdD<-c(ctdD,v)

    }

    ctdC[is.na(ctdC)]<-0

    if(normalized==TRUE){
      ctdC<-ctdC/lenaa
    }
    featureMatrix[n,stC:enC]<-ctdC


    ctdD<-ctdD/lenaa
    featureMatrix[n,stD:enD]<-ctdD

  }

  if(length(label)==numSeqs){
    featureMatrix<-as.data.frame(featureMatrix)
    featureMatrix<-cbind(featureMatrix,label)
  }
  row.names(featureMatrix)<-names(seqs)
  return(featureMatrix)
}
