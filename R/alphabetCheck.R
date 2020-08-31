#' AlphabetCheck
#'
#' This function checks the alphabets in a sequence. If one of the following conditions hold, the sequence will be deleted:
#' 1. A peptide sequence containing non-standard amino acids,
#' 2. A DNA sequence with an alphabet other than A, C, G, or T,
#' 3. An RNA sequence having an alphabet other than A, C, G, or U.
#'
#' @note This function receives a sequence vector and the label of sequences (if any).
#' It deletes sequences (and their labels) containing non-standard alphabets.
#'
#' @param sequences is a string vector. Each element is a peptide, protein, DNA, or RNA sequences.
#'
#' @param alphabet This parameter shows the alphabet of sequences. If it is set to 'aa', it indicates the
#' alphabet of amino acids. When it is 'dna', it shows the nucleotide alphabet and in case it equals 'rna', it represents ribonucleotide alphabet.
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @return 'alphabetCheck' returns a list with two elements. The first element is a vector which contains
#' valid sequences. The second element is a vector which contains the labels of the sequences (if any exists).
#'
#'
#' @export
#'
#' @examples
#'
#' seq<-alphabetCheck(sequences=c("AGDFLIAACNMLKIVYT","ADXVGAJK"),alphabet="aa")
#'

alphabetCheck<-function(sequences,alphabet="aa",label=c()){



  if(length(sequences)==0){
    stop("ERROR: sequence parameter is empty")
  }
  if(length(label)!=0&&length(label)!=length(sequences)){
    stop("ERROR: The lenght of the label vector and the number of sequences do not match!")
  }
  if(alphabet=="rna"){
    alphabet<-c("A","C","G","U")
  }else if (alphabet=="dna"){
    alphabet<-c("A","C","G","T")
  } else if (alphabet=="aa"){
    alphabet<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y" )
  } else {
    stop("ERROR: alphabet shoud be 'dna' or 'rna' or 'aa' ")
  }

  alphabetCheck = sapply(sequences, function(i) all(strsplit(i,
                                                           split = "")[[1]] %in% alphabet))
  flag=0
  if(length(label)==length(sequences)){
    flag=1
    label = label[alphabetCheck]
  } else if(length(label)>0 && length(label)!=length(sequences)){
    stop("ERROR: The number of labels is not equal to the number of sequences!")
  }
  if(is.null(names(sequences))){
    names(sequences)<-as.character(1:length(sequences))
  }
  nonstanSeq<-names(sequences)[!alphabetCheck]

  if(length(nonstanSeq)!=0){
    nonstanSeq<-toString(nonstanSeq)
    warMessage<-paste("The sequences (",nonstanSeq,") were deleted. They contained non-standard alphabets")
    message(warMessage)
  }
  sequences = sequences[alphabetCheck]
  if(length(sequences)==0){
    stop("All sequences contained non-standard alphabets. No sequences remained for analysis :) ")
  }


  if(flag==1){
    names(label)=names(sequences)
  }

  seq_lab<-list("sequences"=sequences,"Lab"=label)

  return(seq_lab)
}
