#' CTD Distribution (CTDD)
#'
#' This function computes the distribution part of \link{CTD}. It calculates fifteen values for each property.
#' For more information, please check the references.
#'
#' @references Dubchak, Inna, et al. "Prediction of protein folding class using global description of amino acid sequence." Proceedings of the National Academy of Sciences 92.19 (1995): 8700-8704.
#'
#' @param seqs is a FASTA file with amino acid sequences. Each sequence starts
#' with a '>' character. Also, seqs could be a string vector. Each element of the vector is a peptide/protein sequence.
#'
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @return This function returns a feature matrix. The number of rows is equal to the number of sequences and the number of columns is 15*7.
#'
#'
#' @export
#'
#' @examples
#'
#' filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
#' CTD_D<-CTDD(seqs=filePrs)
#'

CTDD<-function(seqs,label=c()){

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


  ncharSeqs<-sapply(seqs, nchar)


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
  l=list()
  featureMatrix<-matrix(0,ncol = (numGrp*15),nrow = numSeqs)
  for(n in 1:numSeqs){

    seq=seqs[n]
    numChar=nchar(seq)
    matProp=matrix(0,ncol =numChar,nrow = 7)
    aa=unlist(strsplit(seq,split = ""))

    g1 <- lapply(group1, function(g) which(aa %in% g))

    g2 <- lapply(group2, function(g) which(aa %in% g))

    g3 <- lapply(group3, function(g) which(aa %in% g))

    CTDCvectg1<-lapply(g1, length)
    CTDCvectg2<-lapply(g2, length)
    CTDCvectg3<-lapply(g3, length)
    ctdC<-mapply(c,CTDCvectg1,CTDCvectg2,CTDCvectg3, SIMPLIFY=FALSE)
    ctdC<-unlist(ctdC)
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

    vect<-c()

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
      vect<-c(vect,v)

    }
    numChar
    vect<-vect/numChar
    featureMatrix[n,]<-vect

  }



   nameP1<-paste("P",rep(1:numGrp,each=(3*5)),sep = "")
   nameP2<-paste("G",rep(c(1,2,3),each=5,numGrp),sep="")
   nameP3<-paste(rep(c("First","25%","50%","75%","100%"),(numGrp*3)))
   tempnName<-paste(nameP1,nameP2,nameP3)
   colnames(featureMatrix)<-tempnName


  return(featureMatrix)



}
