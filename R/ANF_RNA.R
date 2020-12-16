#'  Accumulated riboNucleotide Frequency (ANF_RNA)
#'
#' This function replaces ribonucleotides with a four-length vector.
#' The first three elements represent the ribonucleotides and
#' the forth holds the frequency of the ribonucleotide from the beginning of the sequence until the position of the ribonucleotide in the sequence.
#' 'A' will be replaced with c(1, 1, 1, freq), 'C' with c(0, 1, 0, freq),'G' with c(1, 0, 0, freq), and 'U' with c(0, 0, 1, freq).
#'
#' @note This function is provided for sequences with the same lengths.
#' Users can use 'txt' option in outFormat for sequences with different lengths.
#' Warning: If outFormat is set to 'mat' for sequences with different lengths, it returns an error.
#' Also, when output format is 'txt', label information is not shown in the text file.
#' It is noteworthy that 'txt' format is not usable for machine learning purposes if sequences have different sizes. Otherwise 'txt' format
#' is also usable for machine learning purposes.
#'
#' @references Chen, W., Tran, H., Liang, Z. et al. Identification and analysis of the N6-methyladenosine in the Saccharomyces cerevisiae transcriptome. Sci Rep 5, 13859 (2015).
#'
#' @param seqs is a FASTA file containing ribonucleotide sequences. The sequences start
#' with '>'. Also, seqs could be a string vector. Each element of the vector is a ribonucleotide sequence.
#'
#'
#' @param outFormat (output format) can take two values: 'mat'(matrix) and 'txt'. The default value is 'mat'.
#'
#' @param outputFileDist shows the path and name of the 'txt' output file.
#'
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @return The output depends on the outFormat parameter which can be either 'mat' or 'txt'. If outFormat is 'mat', the function returns a feature
#' matrix for sequences with the same length such that the number of columns is (sequence length)*(4)
#' and the number of rows is equal to the number of sequences.
#' If the outFormat is 'txt', the output is written to a tab-delimited file.
#'
#'
#' @export
#'
#' @examples
#'
#' dir = tempdir()
#' fileLNC<-system.file("extdata/Carica_papaya101RNA.txt",package="ftrCOOL")
#' mat<-ANF_RNA(seqs = fileLNC,outFormat="mat")
#'

ANF_RNA<-function(seqs,outFormat="mat",outputFileDist="",label=c()){


  if(length(seqs)==1&&file.exists(seqs)){
    seqs<-fa.read(seqs,alphabet="rna")
    seqs_Lab<-alphabetCheck(seqs,alphabet = "rna",label)

    seqs<-seqs_Lab[[1]]
    label<-seqs_Lab[[2]]
  }
  else if(is.vector(seqs)){
    seqs<-sapply(seqs,toupper)

    seqs_Lab<-alphabetCheck(seqs,alphabet = "rna",label)


    seqs<-seqs_Lab[[1]]
    label<-seqs_Lab[[2]]

  }
  else {
    stop("ERROR: Input sequence is not in the correct format. It should be a FASTA file or a string vector.")
  }

  lenSeqs<-sapply(seqs,nchar)


  nucs<-list("A"=c(1,1,1),"C"=c(0,1,0),"G"=c(1,0,0),"T"=c(0,0,1),"U"=c(0,0,1))
  numSeqs<-length(seqs)

  if(outFormat=="mat"){

    if(length(unique(lenSeqs))>1){
      stop("ERROR: All sequences should have the same length in 'mat' mode. For sequences with different lengths, please use 'txt' for outFormat parameter")
    }


    featureMatrix<-sapply(seqs,function(x){
      chars<-unlist(strsplit(x,""))
      #patA<-`attributes<-`(gregexpr(pattern ="A",x)[[1]],NULL)
      patA<-which(chars=="A")
      lenA<-length(patA)
      tempA<-1:lenA
      valA<-tempA/patA

      #patC<-`attributes<-`(gregexpr(pattern ="C",x)[[1]],NULL)
      patC<-which(chars=="C")
      lenC<-length(patC)
      tempC<-1:lenC
      valC<-tempC/patC

      #patG<-`attributes<-`(gregexpr(pattern ="G",x)[[1]],NULL)
      patG<-which(chars=="G")
      lenG<-length(patG)
      tempG<-1:lenG
      valG<-tempG/patG

      #patT<-`attributes<-`(gregexpr(pattern ="U",x)[[1]],NULL)
      patT<-which(chars=="U")
      lenT<-length(patT)
      tempT<-1:lenT
      valT<-tempT/patT

      vals<-list("A"=valA,"C"=valC,"G"=valG,"U"=valT)
      counter<-list("A"=0,"C"=0,"G"=0,"U"=0)
      len=lenA+lenC+lenG+lenT
      soretdVals<-vector(mode = "numeric",length = len)
      tempVect<-c()
      for(i in 1:len){
        counter[[chars[i]]]=counter[[chars[i]]]+1
        value<-vals[[chars[i]]][counter[[chars[i]]]]
        tempVect<-c(tempVect,nucs[[chars[i]]],value)

      }
      return(tempVect)

    })

    return(t(featureMatrix))

  } else if(outFormat=="txt"){

    nameSeq<-names(seqs)
    featureList<-lapply(seqs,function(x){
      chars<-unlist(strsplit(x,""))
      #patA<-`attributes<-`(gregexpr(pattern ="A",x)[[1]],NULL)
      patA<-which(chars=="A")
      lenA<-length(patA)
      tempA<-1:lenA
      valA<-tempA/patA

      #patC<-`attributes<-`(gregexpr(pattern ="C",x)[[1]],NULL)
      patC<-which(chars=="C")
      lenC<-length(patC)
      tempC<-1:lenC
      valC<-tempC/patC

      #patG<-`attributes<-`(gregexpr(pattern ="G",x)[[1]],NULL)
      patG<-which(chars=="G")
      lenG<-length(patG)
      tempG<-1:lenG
      valG<-tempG/patG

      #patT<-`attributes<-`(gregexpr(pattern ="U",x)[[1]],NULL)
      patT<-which(chars=="U")
      lenT<-length(patT)
      tempT<-1:lenT
      valT<-tempT/patT

      vals<-list("A"=valA,"C"=valC,"G"=valG,"U"=valT)
      counter<-list("A"=0,"C"=0,"G"=0,"U"=0)
      len=lenA+lenC+lenG+lenT
      soretdVals<-vector(mode = "numeric",length = len)
      tempVect<-c()
      for(i in 1:len){
        counter[[chars[i]]]=counter[[chars[i]]]+1
        value<-vals[[chars[i]]][counter[[chars[i]]]]
        tempVect<-c(tempVect,nucs[[chars[i]]],value)

      }
      return(tempVect)

    })
    for(i in 1:numSeqs){
      tem=featureList[[i]]
      temp<-c(nameSeq[i],tem)
      temp<-paste(temp,collapse = "\t")
      write(temp,outputFileDist,append = TRUE)
    }
  }
  else {
    stop("ERROR: outFormat should be 'mat' or 'txt' ")
  }

}
