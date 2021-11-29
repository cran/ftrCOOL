#' Overlapping property features_7bit_T3 (OPF_7bit_T3)
#'
#' This group of functions (OPF Group) categorize amino acids in different groups based on the type.
#' This function includes 7 amino acid properties. OPF_7bit_T3 substitutes each amino acid with a 7-dimensional vector.
#' Each element of the vector shows if that amino acid locates in a special property category or not. '0' means that amino acid is not located in that property group and '1' means it is located.
#' The only difference between OPF_7bit type1, type2, and type3 is in localization of amino acids in the properties groups.
#'
#' @references Wei,L., Zhou,C., Chen,H., Song,J. and Su,R. ACPred-FL: a sequence-based predictor using effective feature representation to improve the prediction of anti-cancer peptides. Bioinformatics (2018).
#'
#' @note This function is provided for sequences with the same lengths.
#' Users can use 'txt' option in outFormat for sequences with different lengths.
#' Warning: If outFormat is set to 'mat' for sequences with different lengths, it returns an error.
#' Also, when output format is 'txt', label information is not shown in the text file.
#' It is noteworthy that 'txt' format is not usable for machine learning purposes if sequences have different sizes. Otherwise 'txt' format
#' is also usable for machine learning purposes.
#'
#' @param seqs is a FASTA file with amino acid sequences. Each sequence starts
#' with a '>' character. Also, seqs could be a string vector. Each element of the vector is a peptide/protein sequence.
#'
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @param outFormat (output format) can take two values: 'mat'(matrix) and 'txt'. The default value is 'mat'.
#'
#' @param outputFileDist shows the path and name of the 'txt' output file.
#'
#'
#' @return The output is different depending on the outFormat parameter ('mat' or 'txt').
#' If outFormat is set to 'mat', it returns a feature matrix for sequences with the same lengths.
#' Number of columns for this feature matrix is equal to (length of the sequences)*7 and number of rows is equal to the number of sequences.
#' If outFormat is 'txt', all binary values will be written to a the output is written to a tab-delimited file. Each line in the file shows the binary format of a sequence.
#'
#'
#' @export
#' @examples
#'
#' ptmSeqsADR<-system.file("extdata/",package="ftrCOOL")
#' ptmSeqsVect<-as.vector(read.csv(paste0(ptmSeqsADR,"/ptmVect101AA.csv"))[,2])
#' mat<-OPF_7bit_T3(seqs = ptmSeqsVect,outFormat="mat")
#'
OPF_7bit_T3<-function(seqs,label=c(),outFormat="mat",outputFileDist=""){

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


  numSeqs<-length(seqs)
  lenSeqs<-sapply(seqs, nchar)

  group<-list("Hydrophobicity"= c('K','R'),
              "Normalized Vander Waals volume"= c('D', 'E', 'K', 'N', 'Q', 'R'),
              "Polarity"=   c('F', 'H', 'K', 'M', 'R', 'W', 'Y'),
              "Polarizibility"=  c('D', 'E', 'H', 'K', 'N', 'Q', 'R'),
              "Charge"= c('F', 'H', 'K', 'M', 'R', 'W', 'Y'),
              "Secondary structures"= c('C', 'F', 'I', 'T', 'V', 'W', 'Y'),
              "Solvent accessibility"=   c('D','E','K','N','R','Q'))


  properties<-c("Hydrophobicity", "Normalized Vander Waals volume",
                "Polarity", "Polarizibility", "Charge", "Secondary structures", "Solvent accessibility")


  if(outFormat=="mat"){
    if(length(unique(lenSeqs))>1){
      stop("ERROR: All sequences should have the same length in 'mat' mode. For sequences with different lengths, please use 'txt' for outFormat parameter")
    }
    featureMatrix<-matrix(0, nrow = numSeqs, ncol = (lenSeqs[1]*7))
    tempN1<-rep(properties,lenSeqs[1])
    tempN2<-rep(1:lenSeqs[1],each=7)
    colnames(featureMatrix)<-paste0("pos",tempN2,"_",tempN1)

    for(n in 1:numSeqs){

      seq=seqs[n]

      aa=unlist(strsplit(seq,split = ""))

      vect<-c()
      for(a in aa)
      {
        g1 <- lapply(group, function(g) which(a %in% g))
        b=lapply(g1, function(x) length(x)>0)
        vect<-c(vect,as.numeric(b))
      }

      featureMatrix[n,]<-vect
    }
    row.names(featureMatrix)<-names(seqs)
    return(featureMatrix)
  }
  else{
    nameSeq<-names(seqs)
    for(n in 1:numSeqs){
      seq<-seqs[n]
      chars<-unlist(strsplit(seq,split = ""))
      vect<-c()
      for(a in aa)
      {
        g1 <- lapply(group, function(g) which(a %in% g))
        b=lapply(g1, function(x) length(x)>0)
        vect<-c(vect,as.numeric(b))
      }
      temp<-c(nameSeq[n],vect)
      temp<-paste(temp,collapse = "\t")
      write(temp,outputFileDist,append = TRUE)

    }
  }

}
