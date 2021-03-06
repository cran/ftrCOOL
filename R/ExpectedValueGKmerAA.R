#' Expected Value for Grouped K-mer Amino Acid(ExpectedValueGKmerAA)
#'
#' This function is introduced by this package for the first time.
#' In this function, amino acids are first grouped into user-defined categories.
#' Later, the expected value of grouped k-mer is computed.
#' Please note that this function differs from Function \link{ExpectedValueKmerAA} which works on individual amino acids.
#'
#'
#' @param seqs is a FASTA file with amino acid sequences. Each sequence starts
#' with a '>' character. Also, seqs could be a string vector. Each element of the vector is a peptide/protein sequence.
#'
#' @param k is an integer. The default value is two.
#'
#' @param Grp is a list of vectors containig amino acids. Each vector represents a category. Users can define a customized amino acid grouping, provided that the sum of all amino acids is 20 and there is no repeated amino acid in the groups.
#' Also, users can choose 'cTriad'(conjointTriad), 'locFus', or 'aromatic'. Each option provides specific information about the type of an amino acid grouping.
#'
#'
#' @param normalized is a logical parameter. When it is FALSE, the return value of the function does not change. Otherwise, the return value is normalized using the length of the sequence.
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#'
#' @return This function returns a feature matrix. The number of rows is equal to the number of sequences and
#' the number of columns is (number of categorizes)^k.
#'
#'
#' @export
#'
#' @examples
#'
#' filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
#' mat1<-ExpectedValueGKmerAA(seqs=filePrs,k=2,Grp="locFus")
#'
#' mat2<-ExpectedValueGKmerAA(seqs=filePrs,k=1,Grp=
#' list(Grp1=c("G","A","V","L","M","I","F","Y","W"),Grp2=c("K","R","H","D","E")
#' ,Grp3=c("S","T","C","P","N","Q")))


ExpectedValueGKmerAA<-function(seqs,k=2,Grp="locFus",normalized=TRUE,label=c())
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
  featureMatrix<-kGAAComposition(seqs,rng=k,Grp=Grp,normalized=FALSE)
  GAAcompos<-kGAAComposition(seqs,rng=1,Grp=Grp,normalized=FALSE)



  for(i in 1:ncol(featureMatrix)){
    chars<-unlist(strsplit(colnames(featureMatrix)[i],split = ""))
    composMat<-GAAcompos[,chars]
    if(is.matrix(composMat)){
      mult<-apply(composMat, 1, prod)
    } else {
      mult<-prod(composMat)
    }

    featureMatrix[,i]<-featureMatrix[,i]/mult
  }
  row.names(featureMatrix)<-names(seqs)

  featureMatrix[is.na(featureMatrix)]=0
  featureMatrix[is.infinite(featureMatrix)]=0

  if(normalized==TRUE){
    seqLen<-sapply(seqs, nchar)
    featureMatrix<-featureMatrix/seqLen
  }

  if(length(label)==numSeqs){
    featureMatrix<-as.data.frame(featureMatrix)
    featureMatrix<-cbind(featureMatrix,label)
  }

  return(featureMatrix)
}
