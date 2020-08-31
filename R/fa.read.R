#' Fasta File Reader
#'
#' This function reads a FASTA file. Each sequence starts with '>' in the file.
#' This is a general function which can be applied to all types of sequences (i.e., protein/peptide, dna, and rna).
#'
#' @param file The address of the FASTA file.
#'
#' @references https://cran.r-project.org/web/packages/rDNAse/index.html
#'
#' @param legacy.mode comments all lines which start with ";".
#'
#' @param seqonly if it is set to true, the function will return sequences with no description.
#'
#' @param alphabet is a vector which contains amino acid, RNA, or DNA alphabets.
#'
#' @return a string vector such that each element is a sequence.
#'
#'
#' @export
#'
#' @examples
#'
#' fileLNC<-system.file("extdata/Athaliana_LNCRNA.fa",package="ftrCOOL")
#' sequenceVectLNC<-fa.read(file=fileLNC,alphabet="dna")
#'
#' filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
#' sequenceVectPRO<-fa.read(file=filePrs,alphabet="aa")


fa.read<-function (file, legacy.mode = TRUE, seqonly = FALSE,alphabet="aa")
{
  if(alphabet=="rna"){
    alphabet<-c("A","C","G","U")
  }else if (alphabet=="dna"){
    alphabet<-c("A","C","G","T")
  } else if (alphabet=="aa"){
    alphabet<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y" )
  } else {
    stop("alphabet shoud be 'dna' or 'rna' or 'aa' ")
  }
  lines = readLines(file)
  if (legacy.mode) {
    comments = grep("^;", lines)
    if (length(comments) > 0) {
      lines = lines[-comments]
    }
  }
  ind = which(substr(lines, 1L, 1L) == ">")
  nseq = length(ind)
  if (nseq == 0)
    stop("There is no line which starts with '>' ")
  start = ind + 1
  end = ind - 1
  end = c(end[-1], length(lines))
  sequences = lapply(seq_len(nseq), function(i) paste(lines[start[i]:end[i]],
                                                      collapse = ""))
  sequences = sapply(sequences, toupper)

  if (seqonly)
    return(sequences)
  nomseq = lapply(seq_len(nseq), function(i) {
    firstword = strsplit(lines[ind[i]], " ")[[1]][1]
    substr(firstword, 2, nchar(firstword))
  })

  names(sequences) = nomseq
  return(sequences)
}



