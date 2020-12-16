#' reverseCompelement (revComp)
#'
#' This function returns the reverse compelement of a dna sequence.
#'
#' @param seq is a dna sequence.
#'
#' @param outputType this parameter can take two values: 'char' or 'str'. If outputType is
#' 'str', the reverse complement sequence of the input sequence is returned as a string.
#' Otherwise, a vector of characters which represent the reverse complement is returned.
#' Default value is 'str'.
#'
#'
#'
#' @return The reverse complement of the input sequence.
#'
#'
#'
#' @export
#'
#' @examples
#'
#' ptmSeqsADR<-system.file("extdata/",package="ftrCOOL")
#' ptmSeqsVect<-as.vector(read.csv(paste0(ptmSeqsADR,"/ptmVect101AA.csv"))[,2])
#' Seq<-ptmSeqsVect[1]
#' revCompSeq<-revComp(seq=Seq,outputType="char")

revComp<-function(seq,outputType="str")
{
  if(outputType!="char"&&outputType!="str"){
    stop("ERROR: outputType should be 'char' or 'str' ")
  }
  temp1<-strsplit(seq,NULL)
  temp2<-rev(unlist(temp1))

  revComp=temp2
  revComp[temp2=="A"]="T"
  revComp[temp2=="T"]="A"
  revComp[temp2=="C"]="G"
  revComp[temp2=="G"]="C"

  if(outputType=="str"){
    revComp=paste(revComp,collapse = "")
  }

  return(revComp)

}
