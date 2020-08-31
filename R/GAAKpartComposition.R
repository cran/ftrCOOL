#' Grouped Amino Acid K Part Composition
#'
#' In this function, amino acids are first grouped into user-defined categories.
#' Later, the composition of the grouped amino acid k part is computed.
#' Please note that this function differs from \link{AAKpartComposition} which works on individual amino acids.
#'
#'
#' @note Warning: The length of all sequences should be greater than k.
#'
#' @param seqs is a FASTA file with amino acid sequences. Each sequence starts
#' with a '>' character. Also, seqs could be a string vector. Each element of the vector is a peptide/protein sequence.
#'
#'
#' @param k is an integer. Each sequence should be divided to k partition(s).
#'
#' @param normalized is a logical parameter. When it is FALSE, the return value of the function does not change. Otherwise, the return value is normalized using the length of the sequence.
#'
#'
#' @param Grp is a list of vectors containig amino acids. Each vector represents a category. Users can define a customized amino acid grouping, provided that the sum of all amino acids is 20 and there is no repeated amino acid in the groups.
#' Also, users can choose 'cTriad'(conjointTriad), 'locFus', or 'aromatic'. Each option provides specific information about the type of an amino acid grouping.
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @return a feature matrix with k*(number of categorizes) number of columns. The number of rows is equal to the number of
#' sequences.
#'
#'
#' @export
#'
#' @examples
#'
#' filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
#' mat1<-GAAKpartComposition(seqs=filePrs,k=5,Grp="aromatic")
#'
#' mat2<-GAAKpartComposition(seqs=filePrs,k=3,normalized=FALSE,Grp=
#' list(Grp1=c("G","A","V","L","M","I","F","Y","W"),Grp2=c("K","R","H","D","E")
#' ,Grp3=c("S","T","C","P","N","Q")))


GAAKpartComposition<- function(seqs,k=5,normalized=TRUE,Grp="locFus",label=c()) {


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

  lenSeqs<-sapply(seqs,nchar)


  if(!all(lenSeqs>=k)){
    deletInd<-which(lenSeqs<k)
    deletedNames<-names(deletInd)
    strNames<-toString(deletedNames)
    lens<-lenSeqs[deletInd]
    strlens<-toString(lens)

    warning(paste("Sequences",strNames,"with lengths",strlens,"were deleted. Their lenghts were smaller than k"))

    if(length(label)==length(lenSeqs)){
      label<-label[-deletInd]
    }
    lenSeqs<-lenSeqs[-deletInd]
    seqs<-seqs[-deletInd]
  }


  if(k<=0)
  {
    stop("k should be greater than 1")
  }


  numSeqs<-length(seqs)


  winSize<-ceiling(lenSeqs/k)

  aa<-vector()
  VectGrp<-c("Grp10"='a',"Grp11"='b',"Grp12"='c',"Grp13"='d',"Grp14"='e',"Grp15"='f',"Grp16"='g',"Grp17"='h',"Grp18"='i',"Grp19"='j',"Grp20"='k')


  for (i in 1:numGrp)
  {
    if(i<10){
      vect<-rep(i,length(Group[[i]]))
      aa<-c(aa,vect)
    } else if(i<=20){
      aa<-c(aa,rep(VectGrp[(i-9)],length(Group[[i]])))
    }
  }

  names(aa)<-grps


  featureMatrix<-matrix(0,ncol = (numGrp*k),nrow = numSeqs)

  tname<-nameKmer(k=1,type = "num",numGrp)
  tname<-rep(tname,k)
  wname<-rep(1:k,each =numGrp)
  coln<-paste(tname,"p",wname,sep = "")

  colnames(featureMatrix)<-coln

  for(n in 1:numSeqs){
    seq<-seqs[n]

    seqChars<-unlist(strsplit(seq,split = ""))


    GrpSeq<-aa[seqChars]





    for(i in 1:k)
    {

      winSeq<-GrpSeq[(((i-1)*winSize[n])+1):(i*winSize[n])]
      twinSeq<-table(winSeq)
      ntwinseq<-names(twinSeq)

      for(j in 1:length(twinSeq))
      {
        colindex<-(i-1)*numGrp+as.numeric(ntwinseq[j])
        featureMatrix[n,colindex]<-twinSeq[j]
      }


    }
  }

  if(normalized==TRUE){
    featureMatrix[,]<-apply(featureMatrix, 2, function(i) i/lenSeqs)
  }

  if(length(label)==numSeqs){
    featureMatrix<-as.data.frame(featureMatrix)
    featureMatrix<-cbind(featureMatrix,label)
  }
  row.names(featureMatrix)<-names(seqs)
  return(featureMatrix)

}



