#' k-Spaced Conjoint Triad (conjointTriadKS)
#'
#' This function calculates the grouped tripeptide composition with conjoint
#' triad grouping type. For each k, it creates a 7^3 feature vector.
#' K is the space between the first and the second amino acids
#' and the second and the third amino acids of the tripeptide.
#'
#'
#' @details A tripeptide with k spaces looks like AA1(ss..s)AA2(ss..s)AA3. AA stands for amino acids and s means space.
#'
#' @note 'upto' is enabled only when rng is a number and not a vector.
#'
#' @param seqs is a FASTA file with amino acid sequences. Each sequence starts
#' with a '>' character. Also, seqs could be a string vector. Each element of the vector is a peptide/protein sequence.
#'
#' @param rng This parameter can be a number or a vector. Each element of the vector shows the number of spaces between the first and the second amino acids
#' and the second and the third amino acids of the tripeptide.
#' For each k in the rng vector, a new vector (whose size is 7^3) is created which contains the frequency of tri-amino acid with k gaps.
#'
#' @param normalized is a logical parameter. When it is FALSE, the return value of the function does not change. Otherwise, the return value is normalized using the length of the sequence.
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @param upto It is a logical parameter. The default value is FALSE. If rng is a number and upto is set to TRUE, rng is converted
#' to a vector with values from 0 to rng.
#'
#' @return This function returns a feature matrix. The number of rows is equal to the number of sequences and
#' the number of columns is (7^3)*(length rng vector).
#'
#'
#' @export
#'
#' @examples
#'
#' filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
#' mat1<-conjointTriadKS(filePrs,rng=2,upto=TRUE,normalized=TRUE)
#'
#' mat2<-conjointTriadKS(filePrs,rng=c(1,3,5))
#'

conjointTriadKS <- function(seqs,rng=3,upto=FALSE,normalized=FALSE,label=c()){


  DictAA<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")

  Group =list(Grp1=c("A","G","V" )
              ,Grp2=c("I","L","F","P"),Grp3=c("Y" ,"M","T","S")
              ,Grp4=c("H","N","Q","W"),Grp5=c("R","K"), Grp6=c("D","E")
              ,Grp7=c("C"))

  #read sequences
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
  numGrp<-length(Group)
  numSeqs=length(seqs)


  aa<-vector()

  for (i in 1:numGrp)
  {
    vect<-rep(i,length(Group[[i]]))
    aa<-c(aa,vect)
  }
  grps<-unlist(Group)
  names(aa)<-grps

  if(upto==TRUE && length(rng)==1){
    l<-length(rng)
    l<-rng[l]
    rng<-0:l
  }


  rng <- sort(rng)
  len<-length(rng)

  featureMatrix <- matrix(0 , ncol = len*((numGrp)^3),nrow = numSeqs)
  Tripep<-nameKmer(k=3,type = "num",numGrp)

  #seprate Tripeptides by a space
  for(i in 1:length(Tripep)){
    Tritemp<-unlist(strsplit(Tripep[i],split = ""))
    Tripep[i]<-paste(Tritemp[1],Tritemp[2],Tritemp[3])
  }


  featName<-vector()
  for(i in 1:len){
    featName<-c(featName,gsub(" ",strrep("s",rng[i]),Tripep))
  }
  featName<-paste0("G(",featName,")")

  colnames(featureMatrix)<-featName

  for(n in 1:numSeqs){
    seq<-seqs[n]
    seqChars<-unlist(strsplit(seq,split = ""))
    GrpSeq<-aa[seqChars]
    lenSeq<-length(GrpSeq)
    for(i in 1:len){
      temp1<-GrpSeq[1:(lenSeq-(2*(rng[i]+1)))]
      temp2<-GrpSeq[((rng[i]+1)+1):(lenSeq-rng[i]-1)]
      temp3<-GrpSeq[((2*(rng[i]+1))+1):(lenSeq)]
      kmers<-paste(temp1,temp2,temp3,sep = "")
      tbkmers<-table(kmers)
      nmtbkmers<-names(tbkmers)

      for(j in 1:length(tbkmers)){
        tmp<-unlist(strsplit(nmtbkmers[j],split = ""))
        index<-(as.numeric(tmp[1])-1)*(numGrp^2)+(as.numeric(tmp[2])-1)*numGrp+as.numeric(tmp[3])
        index<-index+(i-1)*(numGrp^3)
        featureMatrix[n,index]<-tbkmers[j]
      }
    }

  }

  if(normalized==TRUE){
    seqLen<-sapply(seqs, nchar)
    featureMatrix<-featureMatrix/seqLen
  }
  if(length(label)==numSeqs){
    featureMatrix<-as.data.frame(featureMatrix)
    featureMatrix<-cbind(featureMatrix,label)
  }
  row.names(featureMatrix)<-names(seqs)


  return(featureMatrix)
}

