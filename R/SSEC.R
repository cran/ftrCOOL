#' Secondary Structure Elements Composition (SSEC)
#'
#' This function works based on the output of PSIPRED which predicts the secondary structure of the amino acids in a sequence.
#' The output of the PSIPRED is a tab-delimited file which contains the secondary structure in the third column.
#' SSEC returns the frequency of the secondary structures (i.e., Helix, Extended, Coil) of the sequences.
#'
#'
#' @param dirPath Path of the directory which contains all output files of PSIPRED. Each file belongs to a sequence.
#'
#' @return It returns a feature matrix which the number of rows is the number of sequences and
#' the number of columns is 3. The first column shows the number of amino acids which participate
#' in the coil structure. The second column shows the number of amino acids in the extended structure
#' and the last column shows the number of amino acids in the helix structure.
#'
#'
#' @export
#'
#' @examples
#'
#'
#' Predss2dir<-system.file("testForder",package="ftrCOOL")
#' Predss2dir<-paste0(Predss2dir,"/ss2Dir/")
#' mat<-SSEC(Predss2dir)
#'

SSEC<-function(dirPath){

  # source("R/readss2Dir.R")

  ssVectorSimple<-readss2Dir(dirPath)

  matrixPred<-matrix(0,ncol = 3,nrow = length(ssVectorSimple))
  colnames(matrixPred)<-c("C","E","H")
  for(i in 1:length(ssVectorSimple)){
    vect<-ssVectorSimple[[i]]
    temp<-table(vect)
    #temp=temp[c("C","E","H")]
    #temp=temp[!is.na(temp)]
    matrixPred[i,names(temp)]<-temp
  }

  return(matrixPred)
}

