#' Read ss2 predicted Directory (readss2Dir)
#'
#'  This function reads a directory that contains the output files of PSIPRED
#' It gets the directory path as the input and returns a list of vectors. Each vector contains the secondary structure of the amino acids in a peptide/protein sequence.
#'
#' @param dirPath The path of the directory which contains all predss2 files. Each file belongs to a sequence.
#'
#' @return returns a list of vectors with all the predicted secondary structure for each amino acid. The length of the
#' list is the number of files(sequences) and the length of each vector is (length sequence(i))
#'
#'
#' @export
#'
#' @examples
#'
#'
#' PredSS2dir<-system.file("testForder",package="ftrCOOL")
#' PredSS2dir<-paste0(PredSS2dir,"/ss2Dir/")
#' listPredVect<-readss2Dir(PredSS2dir)
#'


readss2Dir<-function(dirPath){
  f=list.files(dirPath)
  fpath<-paste0(dirPath,f)
  listPredVect<-list()
  #ad<-system.file("extdata",package="ftrCOOL")
  ad<-tempdir()
  ad<-paste0(ad,"/tempTable.txt")
  for(i in 1:length(fpath)){
    text<-readLines(fpath[i],skipNul = TRUE)
    st<-which(text=="")
    st<-st+1
    en<-length(text)
    writeLines(text[st:en],ad)
    table<-read.delim(ad,header = FALSE,sep = "")
    table<-as.matrix(table)

    file.remove(ad)
    listPredVect[[i]]<-table[,3]
  }
  return(listPredVect)
}
