#' disorder Simple (DisorderS)
#'
#' This function extracts ordered and disordered amino acids in protein or peptide sequences.
#' The input to the function is provided by VSL2 software. The function
#' represent order amino acids by 'O' and disorder amino acids by 'D'.
#'
#'
#'
#'
#' @note This function is provided for sequences with the same lengths.
#' Users can use 'txt' option in outFormat for sequences with different lengths.
#' Warning: If outFormat is set to 'mat' for sequences with different lengths, it returns an error.
#' Also, when output format is 'txt', label information is not shown in the text file.
#' It is noteworthy that 'txt' format is not usable for machine learning purposes if sequences have different sizes. Otherwise 'txt' format
#' is also usable for machine learning purposes.
#'
#' @param dirPath Path of the directory which contains all output files of VSL2. Each file belongs to a sequence.
#'
#' @param outFormat It can take two values: 'mat' (which stands for matrix) and 'txt'. The default value is 'mat'.
#'
#' @param outputFileDist It shows the path and name of the 'txt' output file.
#'
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
#'
#' PredDisdir<-system.file("testForder",package="ftrCOOL")
#' PredDisdir<-paste0(PredDisdir,"/Disdir/")
#' ad1<-paste0(dir,"/disorderS.txt")
#'
#' DisorderS(PredDisdir, outFormat="txt",outputFileDist=ad1)
#'
#' unlink("dir", recursive = TRUE)



DisorderS<-function(dirPath, outFormat="mat",outputFileDist=""){


  disorderVectorSimple<-readDisDir(dirPath)
  lens<-lapply(disorderVectorSimple, length)

    if(outFormat=="mat"){
      if(length(unique(lens))==1){
        simpleMatrix=matrix("",ncol = lens[[1]],nrow = length(disorderVectorSimple))
        colnames(simpleMatrix)<-paste0("pos",1:lens[[1]])
        for(i in 1:nrow(simpleMatrix)){
          simpleMatrix[i,]=disorderVectorSimple[[i]]
        }
        return(simpleMatrix)
      }
      else{
        stop("Error for sequences with different length outFormat could not set mat, you can set it 'txt'")
      }
    } else {
      for(i in 1:length(disorderVectorSimple)){
        temp<-disorderVectorSimple[[i]]
        tempName<-paste0("sequence",i)
        temp2<-c(tempName,temp)
        vect<-paste(temp2,collapse = "\t")
        write(vect,outputFileDist,append = TRUE)
      }


    }




}
