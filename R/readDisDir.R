#' Read disorder predicted Directory (readDisDir)
#'
#' This function reads a directory that contains the output VSL2 files.
#' It gets the directory path as the input and returns a list of vectors. Each vector includes the disorder/order type for the amino acids of the sequence.
#'
#' @param dirPath the path of a directory which contains all the VSL2 output files.
#'
#' @return a list of vectors with all the predicted disorder/order type for each amino acid. The length of the
#' list is equal to the number of files(sequences) and the length of each vector is the length of the sequence(i).
#'
#'
#' @export
#'
#' @examples
#'
#' PredDisdir<-system.file("testForder",package="ftrCOOL")
#' PredDisdir<-paste0(PredDisdir,"/Disdir/")
#' listPredVect<-readDisDir(PredDisdir)
#'

readDisDir<-function(dirPath){
  f=list.files(dirPath)
  fpath<-paste0(dirPath,f)
  listPredVect<-list()
  #ad<-system.file("extdata",package="ftrCOOL")
  ad=tempdir()
  ad<-paste0(ad,"/tempTable.txt")
  for(i in 1:length(fpath)){
    text<-readLines(fpath[i])
    #text=readLines(predisorder)
    stTemp=which(text=="NO.     RES.    PREDICTION      DISORDER")
    st=stTemp+2
    en=length(text)-1

    writeLines(text[st:en],ad)
    table=read.delim(ad,header = FALSE)
    colnames(table)=c("NO.","RES.","PREDICTION","DISORDER")
    table=as.matrix(table)
    table[,4][table[,4]=="."]="O"
    file.remove(ad)

    listPredVect[[i]]<-table[,4]

  }
  return(listPredVect)
}
