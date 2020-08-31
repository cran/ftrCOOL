#' Splitted Group Amino Acid Composition(SGAAC)
#'
#'
#' In this function, amino acids are first grouped into a user-defined category.
#' Later, the splitted amino Acid composition is computed.
#' Please note that this function differs from \link{SAAC} which works on individual amino acids.
#'
#'
#' @param seqs is a FASTA file with amino acid sequences. Each sequence starts
#' with a '>' character. Also, seqs could be a string vector. Each element of the vector is a peptide/protein sequence.
#'
#'
#' @param k shows which type of amino acid composition applies to the parts.
#' For example, the amino acid composition is applied when k=1 and when k=2, the dipeptide Composition is applied.
#'
#'
#' @param numNterm shows how many amino acids should be considered for N-terminal.
#'
#' @param numCterm shows how many amino acids should be considered for C-terminal.
#'
#'
#' @param Grp is a list of vectors containig amino acids. Each vector represents a category. Users can define a customized amino acid grouping, provided that the sum of all amino acids is 20 and there is no repeated amino acid in the groups.
#' Also, users can choose 'cTriad'(conjointTriad), 'locFus', or 'aromatic'. Each option provides specific information about the type of an amino acid grouping.
#'
#' @param normalized is a logical parameter. When it is FALSE, the return value of the function does not change. Otherwise, the return value is normalized using the length of the sequence.
#'
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#'
#' @return It returns a feature matrix. The number of rows is equal to the number of sequences.
#' The number of columns is 3*((number of groups)^k).
#'
#'
#' @export
#'
#' @examples
#' filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
#' mat<-SGAAC(seqs=filePrs,k=1,numNterm=15,numCterm=15,Grp="aromatic")
#'




SGAAC<- function(seqs,k=1,numNterm=25,numCterm=25,Grp="locFus",normalized=TRUE,label=c()) {

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


  lenSeqs<-sapply(seqs, nchar)
  if(!all(lenSeqs>(numNterm+numCterm))){
    deletInd<-which(lenSeqs<(numNterm+numCterm))
    deletedNames<-names(deletInd)
    strNames<-toString(deletedNames)
    lens<-lenSeqs[deletInd]
    strlens<-toString(lens)

    if(length(deletInd)==length(seqs)){
      stop("ERROR: All sequences were deleted!")
    }

    warning(paste("Sequences",strNames,"with lengths",strlens,"were deleted. The aggregate lenght of Nterminal and Cterminal parts in a sequence should be less than the length of the sequence."))

    if(length(label)==length(lenSeqs)){
      label<-label[-deletInd]
    }
    lenSeqs<-lenSeqs[-deletInd]
    seqs<-seqs[-deletInd]
  }
  numSeqs<-length(seqs)
  ntermSeq<-substr(seqs,1,numNterm)
  ctermSeq<-substr(seqs,(lenSeqs-numCterm+1),lenSeqs)
  midtermSeq<-substr(seqs,(numNterm+1),(lenSeqs-numCterm))
  nTermMat<-kGAAComposition(ntermSeq,rng=k,Grp = Grp,normalized = FALSE)
  cTermMat<-kGAAComposition(ctermSeq,rng=k,Grp = Grp,normalized = FALSE)
  midTermMat<-kGAAComposition(midtermSeq,rng=k,Grp = Grp,normalized = FALSE)
  if(normalized==TRUE){
    nTermMat<-nTermMat/numNterm
    cTermMat<-cTermMat/numCterm
    lenMidle<-(lenSeqs-(numNterm+numCterm))
    midTermMat<-midTermMat/lenMidle
  }
  colnames(nTermMat)<-paste("N_",colnames(nTermMat),sep = "")
  colnames(cTermMat)<-paste("C_",colnames(cTermMat),sep = "")
  colnames(midTermMat)<-paste("M_",colnames(midTermMat),sep = "")

  featureMatrix<-cbind(nTermMat,midTermMat,cTermMat)

  if(length(label)==numSeqs){
    featureMatrix<-as.data.frame(featureMatrix)
    featureMatrix<-cbind(featureMatrix,label)
  }
  row.names(featureMatrix)<-names(seqs)
  return(featureMatrix)


}
