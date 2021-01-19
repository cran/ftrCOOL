#' Read PSSM Directory (readPSSMdir)
#'
#' This function reads a directory that contains the output psi-blast.
#' It gets the directory path as the input and returns a list of vectors. Each vector includes the type for the amino acids of the sequence.
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
#' pssmDir<-system.file("testForder",package="ftrCOOL")
#' pssmDir<-paste0(pssmDir,"/PSSMdir/")
#' listPredVect<-readPSSMdir(pssmDir)
#'

readPSSMdir<-function(dirPath){
  f=list.files(dirPath)
  fpath<-paste0(dirPath,f)
  listPredVect<-list()
  #ad<-system.file("extdata",package="ftrCOOL")
  ad=tempdir()
  ad<-paste0(ad,"/tempTable.txt")
  for(i in 1:length(fpath)){
    text<-readLines(fpath[i])
    #text=readLines(predisorder)
    stTemp=which(text=="           A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V   A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V")
    st=stTemp+1
    revText<-rev(text)
    enTemp=which(revText=="                      K         Lambda")
    en=length(text)-enTemp-1

    writeLines(text[st:en],ad)
    table=read.table(ad,header = FALSE)
    #colnames(table)=c("NO.","RES.","PREDICTION","DISORDER")

    table<-table[,3:22]


    table<-t(table)
    vect=as.vector(table)


    file.remove(ad)
    listPredVect[[i]]<-vect

  }
  return(listPredVect)
}
