#' PseAAC of distance-pairs and reduced alphabet (DistancePair)
#'
#' In this function, first amino acids are grouped into a category which is one of 'cp13', 'cp14', 'cp19', 'cp20'. Users choose one of these terms to categorize amino acids.
#' Then DistancePair function computes frequencies of all grouped residues and also all grouped-paired residues with [0,rng] distance. 'rng'
#' is a parameter which already was set by the user.
#'
#'
#' @param seqs is a FASTA file with amino acid sequences. Each sequence starts
#' with a '>' character. Also, seqs could be a string vector. Each element of the vector is a peptide/protein sequence.
#'
#' @param rng This parameter is a number. It shows maximum number of spaces between amino acid pairs.
#' For each k in the rng vector, a new vector (whose size is (number of categorizes)^2) is created which contains the frequency of pairs with k gaps.
#'
#'
#' @param normalized is a logical parameter. When it is FALSE, the return value of the function does not change. Otherwise, the return value is normalized using the length of the sequence.
#'
#' @param Grp for this parameter users can choose between these items: 'cp13', 'cp14', 'cp19', or 'cp20'.
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @return  This function returns a feature matrix. Row length is equal to the number of sequences and
#' the number of columns is (number of categorizes)+((number of categorizes)^2)*(rng+1).
#'
#'
#' @export
#'
#' @examples
#'
#' filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
#' mat1<-DistancePair(seqs=filePrs,rng=2,Grp="cp14")

DistancePair <- function(seqs,rng=3,normalized=TRUE,Grp="cp14",label=c()){


  upto=TRUE
  DictAA<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")

    if(Grp=="cp13"){
      Group=list(Grp1=c("M","F")
                 ,Grp2=c("I","L"),Grp3=c("V"),Grp4=c("A"),Grp5=c("C")
                 ,Grp6=c("W","Y","Q","H","P"),Grp7=c("G"),Grp8=c("T"),Grp9=c("S"),Grp10=c("N"),Grp11=c("R","K"),Grp12=c("D"),Grp13=c("E"))
    } else if (Grp=="cp14"){

      Group =list(Grp1=c("M","I","V")
                  ,Grp2=c("L"),Grp3=c("F"),Grp4=c("W","Y"),Grp5=c("G")
                  ,Grp6=c("P"),Grp7=c("C"),Grp8=c("A"),Grp9=c("S"),Grp10=c("T"),Grp11=c("N"),Grp12=c("H","R","K","Q"),Grp13=c("E"),Grp14=c("D"))

    } else if (Grp=="cp19"){
      Group =list(Grp1=c("P")
                  ,Grp2=c("G"),Grp3=c("E"),Grp4=c("K"),Grp5=c("R")
                  ,Grp6=c("Q"),Grp7=c("D"),Grp8=c("S"),Grp9=c("N"),Grp10=c("T"),Grp11=c("H"),Grp12=c("C"),Grp13=c("I"),Grp14=c("V"),Grp15=c("W"),Grp16=c("Y","F"),Grp17=c("A"),Grp18=c("L"),Grp19=c("M"))
    }else if (Grp=="cp20"){

      Group =list(Grp1=c("P")
                  ,Grp2=c("G"),Grp3=c("E"),Grp4=c("K"),Grp5=c("R")
                  ,Grp6=c("Q"),Grp7=c("D"),Grp8=c("S"),Grp9=c("N"),Grp10=c("T"),Grp11=c("H"),Grp12=c("C"),Grp13=c("I"),Grp14=c("V"),Grp15=c("W"),Grp16=c("Y"),Grp17=c("A"),Grp18=c("L"),Grp19=c("M"),Grp20=c("F"))

    } else {
      if(!is.list(Grp)){
        stop("ERROR: Grp should be either one of 'cp13', 'cp14','cp19', or 'cp20' or a list containing a valid grouping of amino acids")
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


  #read sequences
  if(length(seqs)==1&&file.exists(seqs)){
    seqs<-fa.read(seqs,alphabet="aa")
    seqs_Lab<-alphabetCheck(seqs,alphabet = "aa",label)

    seqs<-seqs_Lab[[1]]
    label<-seqs_Lab[[2]]

  }else if(is.vector(seqs)){
    seqs<-sapply(seqs,toupper)

    seqs_Lab<-alphabetCheck(seqs,alphabet = "aa",label)

    seqs<-seqs_Lab[[1]]
    label<-seqs_Lab[[2]]

  } else {
    stop("ERROR: Input sequence is not in the correct format. It should be a FASTA file or a string vector.")}


  numSeqs=length(seqs)


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

  if(upto==TRUE && length(rng)==1){
    l<-length(rng)
    l<-rng[l]
    rng<-0:l
  }

  rng <- sort(rng)
  len<-length(rng)

  featureMatrix0<-matrix(0,ncol = numGrp,nrow = numSeqs)
  dipep<-nameKmer(k=2,type = "num",numGrp)

  #seprate dipeptides by a space

  colnames(featureMatrix0)<-nameKmer(k=1,type = "num",numGrp)

  seqChars<-lapply(seqs,function(seq){unlist(strsplit(seq,split = ""))})
  GrpSeqs<-lapply(seqChars,function(seq){aa[seq]})
  featureMatrix<-c()

  for(n in 1:numSeqs){
    seq<-seqs[n]
    sCh<-seqChars[[n]]
    GrpSeq<-GrpSeqs[[n]]
    lenSeq<-length(GrpSeq)
    tabaa<-table(GrpSeq)
    featureMatrix0[n,names(tabaa)]<-tabaa
  }
  featureMatrix<-cbind(featureMatrix,featureMatrix0)
  numcol=numGrp^2

  for(i in 1:len){
    tempMat<-matrix(0,nrow = numSeqs,ncol = numcol)
    colnames(tempMat)<-dipep
    for(n in 1:numSeqs)
    {
      temp1<-GrpSeqs[[n]][1:(lenSeq-rng[i]-1)]
      temp2<-GrpSeqs[[n]][((rng[i]+1)+1):(lenSeq)]
      kmers<-paste(temp1,temp2,sep = "")
      tbkmers<-table(kmers)
      tempMat[n,names(tbkmers)]<-tbkmers
    }
    featureMatrix<-cbind(featureMatrix,tempMat)

  }

  for(i in 1:length(dipep)){
    ditemp<-unlist(strsplit(dipep[i],split = ""))
    dipep[i]<-paste(ditemp[1],ditemp[2])
  }


  featName<-vector()
  for(i in 1:len){
    featName<-c(featName,gsub(" ",strrep("s",rng[i]),dipep))
  }
  featName<-paste0("G(",featName,")")

  colnames(featureMatrix)[(numGrp+1):(ncol(featureMatrix))]<-featName


  if(normalized==TRUE){
    seqLen<-sapply(seqs, nchar)
    featureMatrix<-featureMatrix/(seqLen)
  }
  if(length(label)==numSeqs){
    featureMatrix<-as.data.frame(featureMatrix)
    featureMatrix<-cbind(featureMatrix,label)
  }
  row.names(featureMatrix)<-names(seqs)


  return(featureMatrix)
}




