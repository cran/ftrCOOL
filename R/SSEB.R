#' Secondary Structure Elements Binary (SSEB)
#'
#' This function works based on the output of PSIPRED which predicts the secondary structure of the amino acids in a sequence.
#' The output of the PSIPRED is a tab-delimited file which contains the secondary structure in the third column.
#' SSEB gives a binary number (i.e., '001'='H','010'=E','100'='C') for each amino acid.
#'
#' @note This function is provided for sequences with the same lengths.
#' Users can use 'txt' option in the outFormat parameter for sequences with different lengths.
#' Warning: If the outFormat is set to 'mat' for sequences with different lengths, it returns an error.
#' It is noteworthy that 'txt' format is not usable for machine learning purposes.
#'
#' @details This function converts each amino acid to a 3-bit value, such that 2 bits are 0 and
#' 1 bit is 1. The position of 1 shows the type of the secondary structure of the amino acids in the protein/peptide.
#' In this function, '001' is used to show Helix structure, '010' to show Extended structure and '100' to show coil structure.
#'
#'
#'
#' @param binaryType It can take any of the following values: ('strBin','logicBin','numBin').
#' 'strBin'(String binary): each structure is represented by a string containing 3 characters(0-1). Helix  = "001" , Extended = "010" , coil = "100".
#' 'logicBin'(logical value): Each structure is represented by a vector containing 3 logical entries. Helix  = c(FALSE,FALSE,TRUE) , Extended = c(FALSE,TRUE,FALSE) , Coil = c(TRUE,FALSE,FALSE).
#' 'numBin' (numeric bin): Each structure is represented by a numeric (i.e., integer) vector containing 3 numerals. Helix  = c(0,0,1) , Extended = c(0,1,0) , coil = c(1,0,0).
#'
#'
#'
#' @param dirPath Path of the directory which contains all output files of PSIPRED. Each file belongs to a sequence.
#'
#' @param outFormat It can take two values: 'mat' (which stands for matrix) and 'txt'. The default value is 'mat'.
#'
#' @param outputFileDist It shows the path and name of the 'txt' output file.
#'
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
#' ad<-paste0(dir,"/SSEB.txt")
#'
#' Predss2dir<-system.file("testForder",package="ftrCOOL")
#' Predss2dir<-paste0(Predss2dir,"/ss2Dir/")
#' mat<-SSEB(Predss2dir,binaryType="numBin",outFormat="txt",outputFileDist=ad)
#'
#' unlink("dir", recursive = TRUE)


SSEB<-function(dirPath,binaryType="numBin",outFormat="mat",outputFileDist=""){

  ssVectorSimple<-readss2Dir(dirPath)
  lens<-lapply(ssVectorSimple, length)
  numSeqs<-length(ssVectorSimple)
  dict=list("H"=3,"E"=2,"C"=1)
  if(outFormat=="mat"){
    if(length(unique(lens))>1){
      stop("ERROR all sequences should be in the same lengths in 'mat' mode. Use 'txt' mode for outFormat parameter")
    }else {

    if(binaryType=="strBin"){
      simpleMatrix=matrix("",ncol = lens[[1]],nrow = numSeqs)
      colnames(simpleMatrix)<-paste0("pos",1:lens[[1]])
      for(i in 1:nrow(simpleMatrix)){
        predSS<-ssVectorSimple[[i]]
        predSS[predSS=="H"]="001"
        predSS[predSS=="E"]="010"
        predSS[predSS=="C"]="100"
        simpleMatrix[i,]=predSS
      }

    } else if(binaryType=="logicBin"){

      simpleMatrix<-matrix(FALSE,nrow = numSeqs, ncol = (lens[[1]]*3))
      rng<-(0:(lens[[1]]-1))*3
      for(n in 1:numSeqs){
        predSS<-ssVectorSimple[[n]]
        pos1<-as.numeric(dict[predSS])
        pos1<-rng+pos1
        simpleMatrix[n,pos1]<-TRUE
      }


    } else if(binaryType=="numBin"){
      simpleMatrix<-matrix(0,nrow = numSeqs, ncol = (lens[[1]]*3))
      rng<-(0:(lens[[1]]-1))*3
      for(n in 1:numSeqs){
        predSS<-ssVectorSimple[[n]]
        pos1<-as.numeric(dict[predSS])
        pos1<-rng+pos1
        simpleMatrix[n,pos1]<-1
      }
    } else{
      stop("ERROR choose between these three ('strBin' , 'logicBin' , 'numBin') types of binaryFormat")
    }

      return(simpleMatrix)
    }

  } else {

    for(i in 1:length(ssVectorSimple)){
      predSS<-ssVectorSimple[[i]]
      predSS[predSS=="H"]="001"
      predSS[predSS=="E"]="010"
      predSS[predSS=="C"]="100"

      tempName<-paste0("sequence",i)
      temp<-c(tempName,predSS)
      vect<-paste(temp,collapse = "\t")
      write(vect,outputFileDist,append = TRUE)
    }
  }

}
