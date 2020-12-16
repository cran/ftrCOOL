#'Read Directory of Accessible Solvent accessibility predicted files (readASAdir)
#'
#' This function reads a directory that contains the output files of SPINE-X.
#' It gets the directory path as the input and returns a list of vectors. Each vector includes the ASA predicted value for amino acids of the sequence.
#'
#' @param dirPath path of the directory which contains all the output files of SPINE-X. Each file belongs to a sequence.
#'
#' @return a list of vectors with all the predicted ASA value for each amino acid. The length of the
#' list is the number of files(sequences) and the length of each vector is (length of sequence(i))
#'
#'
#' @export
#'
#' @examples
#'
#'
#'
#' PredASAdir<-system.file("testForder",package="ftrCOOL")
#' PredASAdir<-paste0(PredASAdir,"/ASAdir/")
#' PredVectASA<-readASAdir(PredASAdir)
#'

readASAdir<-function(dirPath){
  f=list.files(dirPath,pattern = )
  fpath<-paste0(dirPath,f)
  listPredVect<-list()
  for(i in 1:length(fpath)){
    table<-read.table(fpath[i],header = FALSE,sep = "",skip = 1)

    table<-table[,11]

    listPredVect[[i]]<-table
  }
  return(listPredVect)
}
