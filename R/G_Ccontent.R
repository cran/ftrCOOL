#' G_C content
#'
#' This function calculates G-C content of each sequence.
#'
#'
#' @param seqs is a FASTA file containing nucleotide sequences. The sequences start
#' with '>'. Also, seqs could be a string vector. Each element of the vector is a nucleotide sequence.
#'
#' @param ORF (Open Reading Frame) is a logical parameter. If it is set to true, ORF region of each sequence is considered instead of the original sequence (i.e., 3-frame).
#'
#' @param reverseORF is a logical parameter. It is enabled only if ORF is true.
#' If reverseORF is true, ORF region will be searched in the sequence and also in the reverse complement of the sequence (i.e., 6-frame).
#'
#'
#' @param normalized is a logical parameter. When it is FALSE, the return value of the function does not change. Otherwise, the return value is normalized using the length of the sequence.
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @return The function returns a feature vector. The length of the vector is equal to the number of sequences.
#' Each entry in the vector contains G-C content of a sequence.
#'
#'
#' @export
#'
#' @examples
#'
#' fileLNC<-system.file("extdata/Athaliana_LNCRNA.fa",package="ftrCOOL")
#' vect<-GCcontent(seqs=fileLNC,ORF=TRUE,reverseORF=FALSE)


GCcontent<-function(seqs,ORF=FALSE,reverseORF=TRUE,normalized=TRUE,label=c()){


  if(length(seqs)==1&&file.exists(seqs)){
    seqs<-fa.read(seqs,alphabet="dna")
    seqs_Lab<-alphabetCheck(seqs,alphabet = "dna",label)

    seqs<-seqs_Lab[[1]]
    label<-seqs_Lab[[2]]
  }
  else if(is.vector(seqs)){
    seqs<-sapply(seqs,toupper)

    seqs_Lab<-alphabetCheck(seqs,alphabet = "dna",label)

    seqs<-seqs_Lab[[1]]
    label<-seqs_Lab[[2]]

  }
  else {
    stop("ERROR: Input sequence is not in the correct format. It should be a FASTA file or a string vector.")
  }

  flag=0
  if(ORF==TRUE){
    if(length(label)==length(seqs)){
      names(label)=names(seqs)
      flag=1
    }
    seqs=maxORF(seqs,reverse=reverseORF)
    if(flag==1)
      label=label[names(seqs)]
  }


  numSeqs<-length(seqs)
  lenSeqs<-sapply(seqs,nchar)
  gcCont<-vector(mode = "numeric",length = numSeqs)
  for(n in 1:numSeqs){
    seq<-seqs[n]
    chars<-unlist(strsplit(seq,""))
    tabChars<-table(chars)
    tabChars<-tabChars[c("A","C","G","T")]
    names(tabChars)<-c("A","C","G","T")
    tabChars[is.na(tabChars)]<-0
    gcCont[n]<-tabChars["G"]+tabChars["C"]

  }
  if(normalized==TRUE){
    gcCont<-gcCont/lenSeqs
  }
  names(gcCont)<-names(seqs)
  if(length(label)==numSeqs){
    gcCont<-as.data.frame(gcCont)
    gcCont<-cbind(gcCont,label)
  }


  return(gcCont)

}


