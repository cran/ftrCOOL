#' riboNucleotide Chemical Property (NCP_RNA)
#'
#' This function replaces ribonucleotides with a three-length vector.
#' The vector represent the ribonucleotides such that
#' 'A' will be replaced with c(1, 1, 1), 'C' with c(0, 1, 0),'G' with c(1, 0, 0), and 'U' with c(0, 0, 1).
#'
#' @note This function is provided for sequences with the same lengths.
#' Users can use 'txt' option in outFormat for sequences with different lengths.
#' Warning: If outFormat is set to 'mat' for sequences with different lengths, it returns an error.
#' Also, when output format is 'txt', label information is not shown in the text file.
#' It is noteworthy that 'txt' format is not usable for machine learning purposes if sequences have different sizes. Otherwise 'txt' format
#' is also usable for machine learning purposes.
#'
#' @references Chen, Zhen, et al. "iLearn: an integrated platform and meta-learner for feature engineering, machine-learning analysis and modeling of DNA, RNA and protein sequence data." Briefings in bioinformatics 21.3 (2020): 1047-1057.
#'
#' @param seqs is a FASTA file containing ribonucleotide sequences. The sequences start
#' with '>'. Also, seqs could be a string vector. Each element of the vector is a ribonucleotide sequence.
#'
#' @param binaryType It can take any of the following values: ('strBin','logicBin','numBin').
#' 'strBin'(String binary): each ribonucleotide is represented by a string containing 4 characters(0-1). A = "0001" , C = "0010" , G = "0100" , T = "1000"
#' 'logicBin'(logical value): Each ribonucleotide is represented by a vector containing 4 logical entries. A = c(F,F,F,T) , ... , T = c(T,F,F,F)
#' 'numBin' (numeric bin): Each ribonucleotide is represented by a numeric (i.e., integer) vector containing 4 numerals. A = c(0,0,0,1) , ... , T = c(1,0,0,0)
#'
#' @param outFormat (output format) can take two values: 'mat'(matrix) and 'txt'. The default value is 'mat'.
#'
#' @param outputFileDist shows the path and name of the 'txt' output file.
#'
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#' @return The output is different depending on the outFormat parameter ('mat' or 'txt').
#' If outFormat is set to 'mat', it returns a feature matrix for sequences with the same lengths.
#' The number of rows is equal to the number of sequences and if binaryType is 'strBin', the number of columns is the length of the sequences.
#' Otherwise, it is equal to (length of the sequences)*3.
#' If outFormat is 'txt', all binary values will be written to a tab-delimited file. Each line in the file shows the binary format of a sequence.
#'
#'
#' @export
#'
#' @examples
#'
#' dir = tempdir()
#' fileLNC<-system.file("extdata/Carica_papaya101RNA.txt",package="ftrCOOL")
#' mat<-NCP_RNA(seqs = fileLNC,binaryType="strBin",outFormat="mat")
#'
#' ad<-paste0(dir,"/NCP.txt")
#' NCP_RNA(seqs = fileLNC,binaryType="numBin",outFormat="txt",outputFileDist=ad)
#' unlink("dir", recursive = TRUE)

NCP_RNA<-function(seqs,binaryType="numBin",outFormat="mat",outputFileDist="",label=c()){


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



  nucs<-list("A"=c(1, 1, 1),"C"=c(0,1,0),"G"=c(1,0,0),"T"=c(0,0,1),"U"=c(0,0,1))
  numSeqs<-length(seqs)

  if(outFormat=="mat"){

    if(length(unique(lenSeqs))>1){
      stop("ERROR: All sequences should have the same length in 'mat' mode. For sequences with different lengths, please use 'txt' for outFormat parameter")
    }

    if(binaryType=="strBin"){

      nucs<-c("A"="111","C"="010","G"="100","T"="001","U"="001")
      featureMatrix<-sapply(seqs,function(x) {
        charList<-unlist(strsplit(x,split = ""))
        cods<-nucs[charList]
        return(cods)
      })
      featureMatrix<-t(featureMatrix)
      colnames(featureMatrix)<-paste("pos",1:lenSeqs[1],sep = "")
      row.names(featureMatrix)<-names(seqs)
    }
    else if(binaryType=="logicBin"){
      nucs<-list("A"=c(TRUE, TRUE, TRUE),"C"=c(FALSE,TRUE,FALSE),"G"=c(TRUE,FALSE,FALSE),"T"=c(FALSE,FALSE,TRUE),"U"=c(FALSE,FALSE,TRUE))
      featureMatrix<-sapply(seqs,function(x) {
        charList<-unlist(strsplit(x,split = ""))
        cods<-nucs[charList]
        cods<-unlist(cods)
        return(cods)
      })
      featureMatrix<-t(featureMatrix)
      #vectCodes<-unlist(cods)
      #featureMatrix<-matrix(vectCodes,byrow = TRUE,ncol = (lenSeqs[1]*3),nrow = numSeqs)
      temp1<-rep(c("P","A","H"),lenSeqs[1])
      temp2<-rep(1:lenSeqs[1],each=3)
      colnames(featureMatrix)<-paste("pos",temp2,"-",temp1,sep = "")
      row.names(featureMatrix)<-names(seqs)
    }
    else if(binaryType=="numBin"){

      featureMatrix<-sapply(seqs,function(x) {
        charList<-unlist(strsplit(x,split = ""))
        cods<-nucs[charList]
        cods<-unlist(cods)
        return(cods)
      })
      featureMatrix<-t(featureMatrix)
      temp1<-rep(c("P","A","H"),lenSeqs[1])
      temp2<-rep(1:lenSeqs[1],each=3)
      colnames(featureMatrix)<-paste("pos",temp2,"-",temp1,sep = "")
      row.names(featureMatrix)<-names(seqs)
    }
    else{
      stop("ERROR! Choose one of 'strBin', 'logicBin', or 'numBin' for binaryFormat")
    }

    return(featureMatrix)

  } else if(outFormat=="txt"){

    nucs<-c("A"="111","C"="010","G"="100","T"="001","U"="001")
    counter<-0
    namesSeqs<-names(seqs)
    codes<-lapply(seqs,function(x) {
      counter<-counter+1
      #print(counter)
      charList<-unlist(strsplit(x,split = ""))
      cods<-nucs[charList]
      namecods<-namesSeqs[counter]
      cods<-unlist(cods)
      cods<-c(namecods,cods)
      temp<-paste(cods,collapse = "\t")
      write(temp,outputFileDist,append = TRUE)
    })
  }
  else {
    stop("ERROR: outFormat should be 'mat' or 'txt' ")
  }

}
