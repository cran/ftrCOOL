#' K Grouped Amino Acid Composition
#'
#' In this function, amino acids are first grouped into user-defined categories.
#' Later, the composition of the k grouped amino acids is computed.
#' Please note that this function differs from \link{kAAComposition} which works on individual amino acids.
#'
#' @details for more details, please refer to \link{kAAComposition}
#'
#' @param seqs is a FASTA file with amino acid sequences. Each sequence starts
#' with a '>' character. Also, seqs could be a string vector. Each element of the vector is a peptide/protein sequence.
#'
#'
#' @param rng This parameter can be a number or a vector. Each entry of the vector holds the value of k in the k-mer composition.
#' For each k in the rng vector, a new vector (whose size is 20^k) is created which contains the frequency of k-mers.
#'
#' @param upto It is a logical parameter. The default value is FALSE. If rng is a number and upto is set to TRUE, rng is converted
#' to a vector with values from 1 to rng.
#'
#'
#' @param normalized is a logical parameter. When it is FALSE, the return value of the function does not change. Otherwise, the return value is normalized using the length of the sequence.
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @param Grp is a list of vectors containig amino acids. Each vector represents a category. Users can define a customized amino acid grouping, provided that the sum of all amino acids is 20 and there is no repeated amino acid in the groups.
#' Also, users can choose 'cTriad'(conjointTriad), 'locFus', or 'aromatic'. Each option provides specific information about the type of an amino acid grouping.
#'
#'
#' @return  This function returns a feature matrix. The number of rows is equal to the number of sequences and
#' the number of columns is ((number of categorizes)^k)*(length of rng vector).
#'
#'
#'
#' @export
#'
#' @examples
#' filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
#' mat1<-CkSGAApair(seqs=filePrs,rng=2,upto=TRUE,Grp="aromatic")
#'
#' mat2<-CkSGAApair(seqs=filePrs,rng=c(1,3,5),Grp=
#' list(Grp1=c("G","A","V","L","M","I","F","Y","W"),Grp2=c("K","R","H","D","E")
#' ,Grp3=c("S","T","C","P","N","Q")))


kGAAComposition<- function(seqs,rng=3,upto=FALSE,normalized=TRUE,Grp="locFus",label=c()) {


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

  }
  else {
    stop("ERROR: Input sequence is not in the correct format. It should be a FASTA file or a string vector.")
  }
  if(upto==TRUE && length(rng)==1){
    l<-length(rng)
    l<-rng[l]
    rng<-0:l
  }
  rng <- sort(rng)
  rng <- unique(rng)
  len<-length(rng)



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

  numSeqs<-length(seqs)

  mergedMatrix<-vector(mode = "numeric")


  for(l in rng){

    featureMatrix<-matrix(0,ncol = ((numGrp)^l),nrow = numSeqs)
    kmername<-nameKmer(l,type="num",num=numGrp)
    colnames(featureMatrix)<-kmername



    for(n in 1:numSeqs){
      seq<-seqs[n]
      seqChars<-unlist(strsplit(seq,split = ""))
      kmers<-""

      GrpSeq<-aa[seqChars]
      lenSeq<-length(GrpSeq)

      #create all kmers occure in the seq
      for (i in 0:(l-1)){
        temp<-GrpSeq[(1+i):(lenSeq-(l-1-i))]
        kmers<-paste(kmers,temp,sep = "")
      }

      # table kmers of the seq
      tabKmers<-table(kmers)

      # a vector with name for each kmer
      tabNames<-names(tabKmers)

      for(i in 1:length(tabKmers))
      {
        temp<-unlist(strsplit(tabNames[i],split = ""))
        num=0
        for(j in 1:l){
          pow<-numGrp^(l-j)
          #num<-num+((as.numeric(temp[j])-1)*pow)
          as.numeric(temp[j])
          num<-(((as.numeric(temp[j]))-1)*pow)+num
        }
        num<-num+1
        featureMatrix[n,num]<-tabKmers[i]
      }

    }

    mergedMatrix<-cbind(mergedMatrix,featureMatrix)

  }

  if(normalized==TRUE){
    seqLen<-sapply(seqs, nchar)
    mergedMatrix<-mergedMatrix/seqLen
  }

  if(length(label)==numSeqs){
    mergedMatrix<-as.data.frame(mergedMatrix)
    mergedMatrix<-cbind(mergedMatrix,label)
  }
  row.names(mergedMatrix)<-names(seqs)
  return(mergedMatrix)

}

