#' Group Dipeptide Deviation from Expected Mean (GrpDDE)
#'
#' This function is introduced by this package for the first time.
#' In this function, amino acids are first grouped into user-defined categories.
#' Later, \link{DDE} is applied to grouped amino acids.
#' Please note that this function differs from \link{DDE} which works on individual amino acids.
#'
#'
#' @param seqs is a FASTA file with amino acid sequences. Each sequence starts
#' with a '>' character. Also, seqs could be a string vector. Each element of the vector is a peptide/protein sequence.
#'
#'
#' @param Grp is a list of vectors containig amino acids. Each vector represents a category. Users can define a customized amino acid grouping, provided that the sum of all amino acids is 20 and there is no repeated amino acid in the groups.
#' Also, users can choose 'cTriad'(conjointTriad), 'locFus', or 'aromatic'. Each option provides specific information about the type of an amino acid grouping.
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @return A feature matrix with (number of categorizes)^2 number of columns. The number of rows is equal to the number of sequences.
#'
#'
#'
#' @export
#'
#' @examples
#' filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
#' mat1<-GrpDDE(seqs=filePrs,Grp="aromatic")
#'
#' mat2<-GrpDDE(seqs=filePrs,Grp=
#' list(Grp1=c("G","A","V","L","M","I","F","Y","W"),Grp2=c("K","R","H","D","E")
#' ,Grp3=c("S","T","C","P","N","Q")))


GrpDDE <-function(seqs,Grp="locFus",label=c())
{

   DictAA<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")

   if(is.list(Grp)){
     Group = Grp
   }
   else{
  if(Grp=="locFus"){
    Group=list(Grp1=c("A","E")
               ,Grp2=c("R","Q","K","H"),Grp3=c("N","D","S","T"),Grp4=c("G"),Grp5=c("P")
               ,Grp6=c("I","L","M","F","V"),Grp7=c("W","Y"),Grp8=c("C"))
  } else if (Grp=="aromatic"){

    Group =list(Grp1=c("G","A","V","L","M","I")
                ,Grp2=c("F","Y","W"),Grp3=c("K","R","H")
                ,Grp4=c("D","E"),Grp5=c("S","T","C","P","N","Q"))

  } else if (Grp=="cTriad"){
    Group =list(Grp1=c("A","G","V" )
                ,Grp2=c("I","L","F","P"),Grp3=c("Y" ,"M","T","S")
                ,Grp4=c("H","N","Q","W"),Grp5=c("R","K"), Grp6=c("D","E")
                ,Grp7=c("C"))
  } else {
    if(!is.list(Grp)){
      stop("ERROR: Grp should be either one of 'locFus', 'aromatic', or 'cTriad' or a list containing a valid grouping of amino acids")
    }

  }
   }


  #Error in inputs (Group members are not unique)
  numGrp<-length(Group)
  grps<-unlist(Group)
  unqGrps<-unique(grps)

  if(!all(grps %in% DictAA)==TRUE){
    stop("ERROR: There is an unknown amino acid in Grp")
  }
  if (length(grps)!=length(unqGrps)){
    stop("ERROR: There is a duplicated amino acid in Grp")
  }
  if(length(grps)!=20)
  {
    stop("ERROR: Total number of amino acids in Grp should be 20 exactly")
  }



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


  #is a matrix with ncol=(lenGrp)^2
  DcGrp <-DcGrp(seqs,grp = Grp)

  # is a vector with length (lenGrp)^2
  TmGrp <-TmGrp(Group)
  names(TmGrp)<-names(DcGrp)

  charSeqs<-sapply(seqs,nchar)

  #is matrix with ncol=(lenGrp)^2
  TvGrp <-TvGrp(TmGrp,charSeqs)


  dde<-t(apply(DcGrp, 1, function(i) i-TmGrp))

  #is a matrix with ncol=(lenGrp)^2
  dde<-dde/sqrt(TvGrp)

  #dde <- (DcGrp-TmGrp)/sqrt(TvGrp)

  if(length(label)==numSeqs){
    dde<-as.data.frame(dde)
    dde<-cbind(dde,label)
  }
  row.names(dde)<-names(seqs)

  return(dde)

}

DcGrp <- function(seqs,grp)
{
  len<-sapply(seqs, nchar)
  len<-len-1
  dipepCompos <- kGAAComposition(seqs,rng=2,normalized = FALSE,Grp = grp)
  dipepCompos <- dipepCompos/len
  return(dipepCompos)
}

TmGrp <- function(grpLst){

  len<-length(grpLst)
  numCodAA <- list("A"=4,"C"=2,"D"=2,"E"=2,"F"=2,"G"=4,"H"=2,"I"=3,"K"=2,"L"=6,"M"=1,"N"=2,"P"=4,"Q"=2,"R"=6,"S"=6,"T"=4,"V"=4,"W"=1,"Y"=2)
  numCodGAA<-vector(mode="numeric",length = len)
  for(i in 1:len){
    sumNum<-0
    for(j in 1: length(grpLst[[i]])){
      sumNum<-sumNum+as.numeric(numCodAA[grpLst[[i]]][j])
    }
    numCodGAA[i]<-sumNum
  }


  Cn<- sum(unlist(numCodAA))
  tm <- matrix(nrow =len ,ncol = len)

  for (i in 1:len){
    for(j in 1:len){
      tm[j,i]<-(numCodGAA[i]*numCodGAA[j])/(Cn*Cn)
    }

  }
  tm<-as.vector(tm)

  return(tm)
}

TvGrp <- function(tm,lens)
{
  x<-tm*(1-tm)
  TvMatrix<-matrix(0,ncol = length(x),nrow = length(lens))
  colnames(TvMatrix)<-colnames(x)
  for(i in 1:length(lens)){
    TvMatrix[i,]<-x/(lens[i]-1)
  }
  return(TvMatrix)

}


