#' Position-Specific Scoring Matrix (PSSM)
#'
#' This functions receives as input PSSM matrices (which are created by PSI-BLAST software)
#' and converts them into feature vectors.
#'
#' @note This function is provided for sequences with the same lengths.
#' Users can use 'txt' option in outFormat for sequences with different lengths.
#' Warning: If outFormat is set to 'mat' for sequences with different lengths, it returns an error.
#' Also, when output format is 'txt', label information is not shown in the text file.
#' It is noteworthy that 'txt' format is not usable for machine learning purposes if sequences have different sizes. Otherwise 'txt' format
#' is also usable for machine learning purposes.
#'
#'
#' @param dirPath Path of the directory which contains all output files of PSI-BLAST. Each file belongs to a sequence.
#'
#' @param outFormat It can take two values: 'mat' (which stands for matrix) and 'txt'. The default value is 'mat'.
#'
#' @param outputFileDist It shows the path and name of the 'txt' output file.
#'
#' @return The output depends on the outFormat parameter which can be either 'mat' or 'txt'. If outFormat is 'mat', the function returns a feature
#' matrix for sequences with the same length such that the number of columns is (sequence length)*(20)
#' and the number of rows is equal to the number of sequences.
#' If the outFormat is 'txt', the output is written to a tab-delimited file.
#'
#'
#' @export
#'
#'
#' @examples
#'
#' dir = tempdir()
#' ad<-paste0(dir,"/pssm.txt")
#'
#' PSSMdir<-system.file("testForder",package="ftrCOOL")
#' PSSMdir<-paste0(PSSMdir,"/PSSMdir/")
#' mat<-PSSM(PSSMdir,outFormat="txt",outputFileDist=ad)
#'
PSSM<-function(dirPath,outFormat="mat",outputFileDist=""){

  VectPSSM<-readPSSMdir(dirPath)
  lens<-lapply(VectPSSM, length)
  aa<-c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")

  if(outFormat=="mat"){
    if(length(unique(lens))==1){
      simpleMatrix=matrix("",ncol = lens[[1]],nrow = length(VectPSSM))
      tempname<-1:(lens[[1]]/20)
      tempname<-rep(tempname,each=20)
      colnames(simpleMatrix)<-paste0(rep(aa,(lens[[1]]/20)),"_pos",tempname)
      for(i in 1:nrow(simpleMatrix)){
        simpleMatrix[i,]=VectPSSM[[i]]
      }
      return(simpleMatrix)
    }
    else{
      stop("ERROR all sequences should be in the same lengths in 'mat' mode. Use 'txt' mode for outFormat parameter")
    }
  } else {
    for(i in 1:length(VectPSSM)){
      vectSim<-(VectPSSM[[i]])
      tempName<-paste0("sequence",i)
      temp<-c(tempName,vectSim)
      vect<-paste(temp,collapse = "\t")
      write(vect,outputFileDist,append = TRUE)
    }

  }
}

