#' Read Directory of Torsion predicted files (readTorsionDir)
#'
#' This function reads a directory that contains the output files of SPINE-X.
#' It gets the directory path as the input and returns a list of vectors. Each vector includes the phi and psi angle
#' of the amino acids of the sequence.
#'
#'
#' @param dirPath The path of the directory which contains all output files of SPINE-X. Each file belongs to a sequence.
#'
#'
#' @return returns a list of vectors with all the predicted phi and psi angles for each amino acid. The length of the
#' list is the number of files(sequences) and the length of each vector is (2(phi-psi)*length sequence(i)).
#'
#'
#' @export
#'
#' @examples
#'
#'
#' PredTorsioNdir<-system.file("testForder",package="ftrCOOL")
#' PredTorsioNdir<-paste0(PredTorsioNdir,"/TorsioNdir/")
#' PredVectASA<-readTorsionDir(PredTorsioNdir)
#'


readTorsionDir<-function(dirPath){
  f=list.files(dirPath)
  fpath<-paste0(dirPath,f)
  listPredVect<-list()
  for(i in 1:length(fpath)){
    table<-read.table(fpath[i],header = FALSE,sep = "",skip = 1)


    table<-table[,c(4,5)]


    table<-t(table)
    vect<-as.vector(table)
    listPredVect[[i]]<-vect
  }
  return(listPredVect)
}
