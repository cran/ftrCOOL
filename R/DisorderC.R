#' disorder Content (DisorderC)
#'
#' This function extracts ordered and disordered amino acids in protein or peptide sequences. The input to the function is provided by VSL2 software.
#' Also, the function returns number of order and disorder amino acids in the sequence.
#'
#'
#'
#' @param dirPath Path of the directory which contains all output files of VSL2. Each file belongs to a sequence.
#'
#' @return The output is a feature matrix with 2 columns. The number of rows is equal to the number of sequences.
#'
#'
#' @export
#'
#' @examples
#'
#' dir = tempdir()
#' PredDisdir<-system.file("testForder",package="ftrCOOL")
#' PredDisdir<-paste0(PredDisdir,"/Disdir/")
#'
#' mat<-DisorderC(PredDisdir)
#'



DisorderC<-function(dirPath){


  disorderVectorSimple<-readDisDir(dirPath)

  matrixPred<-matrix(0,ncol = 2,nrow = length(disorderVectorSimple))
  colnames(matrixPred)<-c("D","O")
  for(i in 1:length(disorderVectorSimple)){
    vect<-disorderVectorSimple[[i]]
    temp<-table(vect)
    #temp=temp[c("D","E")]
    #temp=temp[!is.na(temp)]
    matrixPred[i,names(temp)]<-temp
    }
    return(matrixPred)

}
