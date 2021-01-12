#' Accessible Solvent Accessibility (ASA)
#'
#' ASA represents an amino acid by a numeric value.
#' This function extracts the ASA from the output of SPINE-X software which predicts ASA for each amino acid in a peptide or protein sequence.
#' The output of SPINE-X is a tab-delimited file. ASAs are in the 11th column of the file.
#'
#' @note This function is provided for sequences with the same lengths.
#' Users can use 'txt' option in outFormat for sequences with different lengths.
#' Warning: If outFormat is set to 'mat' for sequences with different lengths, it returns an error.
#' Also, when output format is 'txt', label information is not shown in the text file.
#' It is noteworthy that 'txt' format is not usable for machine learning purposes if sequences have different sizes. Otherwise 'txt' format
#' is also usable for machine learning purposes.
#'
#'
#' @param dirPath Path of the directory which contains all output files of SPINE-X. Each file belongs to a sequence.
#'
#' @param outFormat It can take two values: 'mat' (which stands for matrix) and 'txt'. The default value is 'mat'.
#'
#' @param outputFileDist It shows the path and name of the 'txt' output file.
#'
#' @return The output depends on the outFormat which can be either 'mat' or 'txt'. If outFormat is 'mat', the function returns a feature
#' matrix for sequences with the same lengths such that the number of columns is equal to the length of the sequences
#' and the number of rows is equal to the number of sequences.
#' If the outFormat is 'txt', the output is written to a tab-delimited file.
#'
#'
#' @export
#'
#' @import utils
#'
#' @import stats
#'
#' @examples
#'
#' dir = tempdir()
#' ad<-paste0(dir,"/asa.txt")
#'
#'
#' PredASAdir<-system.file("testForder",package="ftrCOOL")
#' PredASAdir<-paste0(PredASAdir,"/ASAdir/")
#' ASA(PredASAdir,outFormat="txt",outputFileDist=ad)
#'
#' unlink("dir", recursive = TRUE)


ASA<-function(dirPath,outFormat="mat",outputFileDist=""){

  VectorSimple<-readASAdir(dirPath)
  lens<-lapply(VectorSimple, length)

  if(outFormat=="mat"){
    if(length(unique(lens))==1){
      simpleMatrix=matrix("",ncol = lens[[1]],nrow = length(VectorSimple))
      colnames(simpleMatrix)<-paste0("pos",1:lens[[1]])
      for(i in 1:nrow(simpleMatrix)){
        simpleMatrix[i,]=VectorSimple[[i]]
      }
      return(simpleMatrix)
    }
    else{
      stop("ERROR all sequences should be in the same lengths in 'mat' mode. Use 'txt' mode for outFormat parameter")
    }
  } else {
    for(i in 1:length(VectorSimple)){
      vectSim<-(VectorSimple[[i]])
      tempName<-paste0("sequence",i)
      temp<-c(tempName,vectSim)
      vect<-paste(temp,collapse = "\t")
      write(vect,outputFileDist,append = TRUE)
    }

  }

}
