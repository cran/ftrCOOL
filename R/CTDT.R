#' CTD Transition (CTDT)
#'
#' This function computes the transition part of \link{CTD}. Thirteen properties are defined in this function. Each
#' property categorizes the amino acids of a sequence into three groups. For each property, the grouped amino acid transition (i.e., transitions 1-2, 1-3, and 2-3)
#' is calculated. For more information, please check the references.
#'
#' @references Dubchak, Inna, et al. "Prediction of protein folding class using global description of amino acid sequence." Proceedings of the National Academy of Sciences 92.19 (1995): 8700-8704.
#'
#'
#' @param seqs is a FASTA file with amino acid sequences. Each sequence starts
#' with a '>' character. Also, seqs could be a string vector. Each element of the vector is a peptide/protein sequence.
#'
#' @param normalized is a logical parameter. When it is FALSE, the return value of the function does not change. Otherwise, the return value is normalized using the length of the sequence.
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @return This function returns a feature matrix.
#' The number of rows is equal to the number of sequences and the number of columns is 3*7,
#' where three is the number of transition types (i.e., 1-2, 1-3, and 2-3) and thirteen is the number of properties.
#'
#'
#' @export
#'
#' @examples
#'
#' filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
#' CTD_T<-CTDT(seqs=filePrs,normalized=FALSE)
#'

CTDT<-function(seqs,normalized=FALSE,label=c()){



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


  numGrp=7

  numSeqs<-length(seqs)

  ctdt<-matrix(0,ncol = 3*7,nrow = numSeqs)
  colnames(ctdt)<-paste0("G",rep(1:numGrp,each=3),rep(c("1,2","1,3","2,3"),7))


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


  for(n in 1:numSeqs){

    seq=seqs[n]

    aa=unlist(strsplit(seq,split = ""))

    g1 <- lapply(group1, function(g) which(aa %in% g))

    g2 <- lapply(group2, function(g) which(aa %in% g))

    g3 <- lapply(group3, function(g) which(aa %in% g))



    G <- vector("list", numGrp)
    lenaa=length(aa)
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
    Vmtx<-as.vector(mtx)
    #Vmtx<-Vmtx/lenaa
    ctdt[n,]<-Vmtx

  }

  FVnameP1<-paste("P",rep(1:numGrp,each=3),sep = "")
  FVnameP2<-paste("T(",rep(c("1,2","1,3","2,3"),numGrp),")",sep="")
  FVname<-paste(FVnameP1,FVnameP2)
  colnames(ctdt)<-FVname

  if(length(label)==numSeqs){
    ctdt<-as.data.frame(ctdt)
    ctdt<-cbind(ctdt,label)
  }
  if(normalized==TRUE){
    seqLen<-sapply(seqs, nchar)
    ctdt<-ctdt/seqLen
  }

  row.names(ctdt)<-names(seqs)


  return(ctdt)



}
