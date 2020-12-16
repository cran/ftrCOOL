#' Secondary Structure Elements Simple (SSES)
#'
#' This function works based on the output of PSIPRED which predicts the secondary structure of the amino acids in a sequence.
#' The output of the PSIPRED is a tab-delimited file which contains the secondary structure in the third column.
#' The function represent amino acids in the helix structure by 'H', amino acids in the extended structure by 'E', and amino acids in the coil structure by 'C'.
#'
#' @note This function is provided for the sequences with the same lengths. However,
#' the users can use 'txt' option in the outFormat parameter for sequences with different lengths.
#' Warning: If the outFormat is set to 'mat' for sequences with different lengths, it returns an error.
#' Also, when the output format is 'txt', the label information is not displayed in the text file.
#' It is noteworthy that, 'txt' format is not usable for machine learning purposes.
#'
#' @param dirPath Path of the directory which contains all output files of PSIPRED. Each file belongs to a sequence.
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
#' @examples
#'
#' dir = tempdir()
#' ad<-paste0(dir,"/simpleSSE.txt")
#'
#' Predss2dir<-system.file("testForder",package="ftrCOOL")
#' Predss2dir<-paste0(Predss2dir,"/ss2Dir/")
#' mat<-SSES(Predss2dir,outFormat="txt",outputFileDist=ad)
#'

SSES<-function(dirPath,outFormat="mat",outputFileDist=""){


  ssVectorSimple<-readss2Dir(dirPath)
  lens<-lapply(ssVectorSimple, length)


  if(outFormat=="mat"){
    if(length(unique(lens))==1){
      simpleMatrix=matrix("",ncol = lens[[1]],nrow = length(ssVectorSimple))
      colnames(simpleMatrix)<-paste0("pos",1:lens[[1]])
      for(i in 1:nrow(simpleMatrix)){
        simpleMatrix[i,]=ssVectorSimple[[i]]
      }
      return(simpleMatrix)
    }
    else{
      stop("ERROR all sequences should be in the same lengths in 'mat' mode. Use 'txt' mode for outFormat parameter")
    }
  } else {
    for(i in 1:length(ssVectorSimple)){

      predSS<-ssVectorSimple[[i]]
      tempName<-paste0("sequence",i)
      temp<-c(tempName,predSS)
      vect<-paste(temp,collapse = "\t")
      write(vect,outputFileDist,append = TRUE)

    }
  }

}
