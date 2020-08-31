#' nonStandard amino acid sequence
#'
#' This function returns sequences which contain at least one non-standard alphabet.
#'
#'
#' @param file The address of fasta file which contains all the sequences.
#'
#' @param legacy.mode It comments all lines starting with ";"
#'
#' @param seqonly If it is set to true, the function returns sequences with no description.
#'
#' @param alphabet It is a vector which contains the amino acid, RNA, or DNA alphabets.
#'
#' @return This function returns a string vector. Each element of the vector is a sequence which contains
#' at least one non-standard alphabet.
#'
#'
#' @export
#'
#' @examples
#'
#' filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
#' nonStandardPrSeq<-nonStandardSeq(file = filePrs,alphabet="aa")
#'
#' fileLNC<-system.file("extdata/Athaliana_LNCRNA.fa",package="ftrCOOL")
#' nonStandardNUCSeq<-nonStandardSeq(file = filePrs, alphabet="dna")


nonStandardSeq<-function (file, legacy.mode = TRUE, seqonly = FALSE,alphabet="aa")
{
  if(alphabet=="rna"){
    alphabet<-c("A","C","G","U")
  }else if (alphabet=="dna"){
    alphabet<-c("A","C","G","T")
  } else if (alphabet=="aa"){
    alphabet<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y" )
  } else {
    stop("ERROR: alphabet parameter shoud set 'dna' or 'rna' or 'aa'")
  }
  if(length(file)==1&&file.exists(file)){
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
    stop("ERROR: There is no line which starts with '>' ")
  start = ind + 1
  end = ind - 1
  end = c(end[-1], length(lines))
  sequences = lapply(seq_len(nseq), function(i) paste(lines[start[i]:end[i]],
                                                      collapse = ""))
  sequences = lapply(sequences, toupper)
  alphabetCheck = sapply(sequences, function(i) all(strsplit(i,
                                                             split = "")[[1]] %in% alphabet))
  nonstandardindex<-which(alphabetCheck==FALSE,arr.ind = TRUE)

  sequences = sequences[nonstandardindex]
  if (seqonly)
    return(sequences)
  nomseq = lapply(seq_len(nseq), function(i) {
    firstword = strsplit(lines[ind[i]], " ")[[1]][1]
    substr(firstword, 2, nchar(firstword))
  })
  names(sequences) = nomseq[nonstandardindex]
  return(sequences)}


  else if(is.vector(file)){


    sequences = lapply(file, toupper)
    alphabetCheck = sapply(sequences, function(i) all(strsplit(i,
                                                               split = "")[[1]] %in% alphabet))
    nonstandardindex<-which(alphabetCheck==FALSE,arr.ind = TRUE)
    if(!is.null(names(sequences))){
      tempName<-names(sequences)[nonstandardindex]
    } else {
      names(sequences)<-as.character(1:length(sequences))
      tempName<-as.character(nonstandardindex)
    }

    sequences = sequences[nonstandardindex]

    names(sequences)<-tempName
    sequences<-unlist(sequences)
    return(sequences)


  } else {
    stop("ERROR: Sequences should be a vector or in fasta file")
  }
}



