#' Conjoint Triad
#'
#' This function calculates the grouped tripeptide composition with the conjoint triad grouping type.
#'
#'
#' @param seqs is a FASTA file with amino acid sequences. Each sequence starts
#' with a '>' character. Also, seqs could be a string vector. Each element of the vector is a peptide/protein sequence.
#'
#' @param normalized is a logical parameter. When it is FALSE, the return value of the function does not change. Otherwise, the return value is normalized using the length of the sequence.
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @return This function returns a feature matrix. The number of rows equals to the number of sequences and
#' the number of columns is 7^3.
#'
#'
#' @export
#'
#' @examples
#'
#' filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
#' mat1<-conjointTriad(seqs=filePrs)
#'

conjointTriad<- function(seqs,normalized=TRUE,label=c()){


cTriad<-kGAAComposition(seqs,rng=3,normalized = normalized,upto = FALSE,Grp = "cTriad",label = label)
  return(cTriad)
}
