#' Enhanced Grouped Amino Acid Composition
#'
#' In this function, amino acids are first grouped into user-defined categories. Then, enhanced grouped amino acid composition is computed.
#' For details about the enhanced feature, please refer to function \link{EAAComposition}.
#' Please note that this function differs from function \link{EAAComposition} which works on individual amino acids.
#'
#'
#' @param seqs is a FASTA file with amino acid sequences. Each sequence starts
#' with a '>' character. Also, seqs could be a string vector. Each element of the vector is a peptide/protein sequence.
#'
#' @param winSize shows the size of sliding window. It is a numeric value.
#'
#' @param overLap This parameter shows how the window
#' moves on the sequence. If the overlap is set to TRUE, the next window would have distance 1 with
#' the previous window. Otherwise, the next window will start from the next amino acid after the previous window.
#' There is no overlap between the next and previous windows.
#'
#' @param Grp is a list of vectors containig amino acids. Each vector represents a category. Users can define a customized amino acid grouping, provided that the sum of all amino acids is 20 and there is no repeated amino acid in the groups.
#' Also, users can choose 'cTriad'(conjointTriad), 'locFus', or 'aromatic'. Each option provides specific information about the type of an amino acid grouping.
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @param outFormat (output format) can take two values: 'mat'(matrix) and 'txt'. The default value is 'mat'.
#'
#' @param outputFileDist shows the path and name of the 'txt' output file.
#'
#' @return The output depends on the outFormat parameter which can be either 'mat' or 'txt'. If outFormat is 'mat', the function returns a feature
#' matrix for sequences with the same length such that the number of columns is ((number of categorizes) * (number of windows))
#' and the number of rows is equal to the number of sequences. It is usable for machine learning purposes.
#' If the outFormat is 'txt', the output is written to a tab-delimited file.
#'
#'
#'
#' @export
#'
#' @examples
#'
#' dir = tempdir()
#' ptmSeqsADR<-system.file("extdata/",package="ftrCOOL")
#' ptmSeqsVect<-as.vector(read.csv(paste0(ptmSeqsADR,"/ptmVect101AA.csv"))[,2])
#' mat1<-EGAAComposition(seqs = ptmSeqsVect,winSize=20,overLap=FALSE,Grp="locFus")
#'
#' mat2<-EGAAComposition(seqs = ptmSeqsVect,winSize=30,overLap=FALSE,Grp=
#' list(Grp1=c("G","A","V","L","M","I","F","Y","W"),Grp2=c("K","R","H","D","E")
#' ,Grp3=c("S","T","C","P","N","Q")),outFormat="mat")
#'
#' ad<-paste0(dir,"/EGrpaaCompos.txt")
#' filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
#' EGAAComposition(seqs = filePrs,winSize=20,Grp="cTriad",outFormat="txt"
#' ,outputFileDist=ad)
#'


EGAAComposition <- function(seqs,winSize=50,overLap=TRUE,Grp="locFus",label=c(),outFormat="mat",outputFileDist="")
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

  lenSeqs<-sapply(seqs,nchar)
  if(!all(lenSeqs>=winSize)){
    deletInd<-which(lenSeqs<winSize)
    deletedNames<-names(deletInd)
    strNames<-toString(deletedNames)
    lens<-lenSeqs[deletInd]
    strlens<-toString(lens)

    warning(paste("Sequences",strNames,"with lengths",strlens,"were deleted. Their lenghts were smaller than the window size"))

    if(length(label)==length(lenSeqs)){
      label<-label[-deletInd]
    }
    lenSeqs<-lenSeqs[-deletInd]
    seqs<-seqs[-deletInd]
  }


  if(overLap==FALSE){
    st<-0:(ceiling(lenSeqs[1]/winSize)-1)
    st<-(st*winSize)+1
    en<-st+winSize-1
    len=ceiling(lenSeqs[1]/winSize)
  }
  else{
    st<-seq(1,(lenSeqs[1]-winSize+1),1)
    en<-seq(winSize,lenSeqs[1],1)
    len<-lenSeqs[1]-winSize+1
  }
  numSeqs<-length(seqs)

  if(outFormat=="mat"){
    if(length(unique(lenSeqs))>1){
      stop("ERROR: All sequences should have the same length in 'mat' mode. For sequences with different lengths, please use 'txt' for outFormat parameter")
    }



    featureMatrix<-matrix(0,ncol = (numGrp*len),nrow = numSeqs)

    tname<-nameKmer(k=1,type = "num",num = numGrp)
    tname<-rep(tname,len)
    wname<-rep(0:(len-1),each =numGrp)

    coln<-paste("G",tname,"-w",wname,sep = "")

    colnames(featureMatrix)<-coln


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


    for(n in 1:numSeqs){

      seq<-seqs[n]
      subSeqs<-substring(seq,st,en)
      temp<-sapply(subSeqs, function(i) kGAAComposition(i,rng = 1,upto = FALSE,normalized = FALSE,Grp=Grp))
      temp<-as.vector(temp)
      featureMatrix[n,]<-temp
    }


    if(length(label)==numSeqs){
      featureMatrix<-as.data.frame(featureMatrix)
      featureMatrix<-cbind(featureMatrix,label)
    }
    row.names(featureMatrix)<-names(seqs)
    return(featureMatrix)


  } else if(outFormat=="txt"){


    aa<-vector()
    nameSeq<-names(seqs)
    for (i in 1:numGrp)
    {
      vect<-rep(i,length(Group[[i]]))
      aa<-c(aa,vect)
    }
    names(aa)<-grps


    for(n in 1:numSeqs){

      seq<-seqs[n]
      subSeqs<-substring(seq,st,en)
      temp<-sapply(subSeqs, function(i) kGAAComposition(i,rng = 1,upto = FALSE,normalized = FALSE,Grp=Grp))
      temp<-unlist(temp)
      temp<-c(nameSeq[n],temp)
      temp2<-paste(temp,collapse = "\t")
      write(temp2,outputFileDist,append = TRUE)
    }


  } else{
    stop("ERROR: outFormat should be 'mat' or 'txt' ")
  }


}



