#' Effective Number of Codon (EffectiveNumberCodon)
#'
#' This function calculates the effective number of codon for each sequence.
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
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#'
#' @return The function returns a feature vector. The length of the vector is equal to the number of sequences.
#' Each entry in the vector contains the effective number of codon.
#'
#'
#' @export
#'
#' @examples
#'
#'
#' fileLNC<-system.file("extdata/Athaliana_LNCRNA.fa",package="ftrCOOL")
#' vect<-EffectiveNumberCodon(seqs=fileLNC,ORF=TRUE,reverseORF=FALSE)
#'


EffectiveNumberCodon<-function(seqs,ORF=FALSE,reverseORF=TRUE,label=c())
{

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
  codonList<-list("A"=c("GCT","GCC","GCA","GCG"),"R"=c("AGA","AGG","CGT","CGC","CGA","CGG"),
                  "N"=c("AAT","AAC"),"D"=c("GAT","GAC"),"C"=c("TGT","TGC"),"Q"=c("CAA","CAG"),
                  "E"=c("GAA","GAG"),"G"=c("GGT","GGC","GGA","GGG"),"H"=c("CAT","CAC"),
                  "I"=c("ATT","ATC","ATA"),"L"=c("TTA","TTG","CTT","CTC","CTA","CTG"),"K"=c("AAA","AAG"),
                  "M"=c("ATG"),"F"=c("TTT","TTC"),"P"=c("CCT","CCC","CCA","CCG"),
                  "S"=c("AGT","AGC","TCT","TCC","TCA","TCG"),"T"=c("ACT","ACC","ACA","ACG"),"W"=c("TGG"),
                  "Y"=c("TAT","TAC"),"V"=c("GTT","GTC","GTA","GTG"),"STOPcod"=c("TAA","TAG","TGA"))
  index<-vector()
  for(i in 1:length(codonList)){
    index<-c(index,length(codonList[[i]]))
  }
  index<-cumsum(index)
  codonVector<-unlist(codonList)
  lenCodList<-length(codonList)

  numSeqs<-length(seqs)
  lenSeqs<-sapply(seqs,nchar)


  featureMatrix<-matrix(0,nrow = numSeqs,ncol = 4^3)
  colnames(featureMatrix)<-codonVector
  effectiveNCVect<-vector(mode = "numeric",length = numSeqs)


  for(n in 1:numSeqs){
    seq<-seqs[n]
    triplet<-substring(seq,seq(1,(lenSeqs[n]-2),3),seq(3,lenSeqs[n],3))
    triplet<-triplet[1:(floor(lenSeqs[n]/3))]
    tableContent<-table(triplet)
    featureMatrix[n,names(tableContent)]<-tableContent
    Nc<-0
    st<-1
    for(j in 1:lenCodList){
      if(j>1){
        st<-index[(j-1)]+1
      }
      end<-index[j]
      tempVect<-featureMatrix[n,st:end]

      m=sum(tempVect)
      if(m>0){
        p=tempVect/m
        p2<-p^2
        fscore<-((m*(sum(p2)))-1)/(m-1)
        if(is.na(fscore)){
          fscore=1
          Nc=Nc+1
        }
        else if(fscore>0){
          Nc=Nc+1/fscore
        }
      }

    }

    effectiveNCVect[n]<-Nc
  }

  names(effectiveNCVect)<-names(seqs)
  if(length(label)==numSeqs){

    effectiveNCVect<-cbind(effectiveNCVect,label)
  }

  return(effectiveNCVect)

}
