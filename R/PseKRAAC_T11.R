#' Pseudo K_tuple Reduced Amino Acid Composition Type-11 (PseKRAAC_T11)
#'
#' There are 16 types of PseKRAAC function.
#' In the functions, a (user-selected) grouping of the amino acids might be used to reduce the amino acid alphabet.
#' Also, the functions have a type parameter.
#' The parameter determines the protein sequence analyses which can be either gap or lambda-correlation.
#' PseKRAAC_type11(PseKRAAC_T11) contains Grp 2-20.
#'
#'
#' @references Zuo, Yongchun, et al. "PseKRAAC: a flexible web server for generating pseudo K-tuple reduced amino acids composition." Bioinformatics 33.1 (2017): 122-124.
#'
#' @param seqs is a FASTA file with amino acid sequences. Each sequence starts
#' with a '>' character. Also, seqs could be a string vector. Each element of the vector is a peptide/protein sequence.
#'
#' @param type This parameter has two valid value "lambda" and "gap". "lambda"
#' calls lambda_model function and "gap" calls gap_model function.
#'
#' @param Grp is a numeric value. It shows the id of an amino acid group.
#' Please find the available groups in the detail section.
#'
#'
#' @param GapOrLambdaValue is an integer.
#' If type is gap, this value shows number of gaps between two k-mers.
#' If type is lambda, the value of GapOrLambdaValue shows the number of gaps between each two amino acids of k-mers.
#'
#' @param k This parameter keeps the value of k in k-mer.
#'
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#'
#' @return This function returns a feature matrix. The number of rows is equal to the number of sequences and
#' the number of columns is (Grp)^k.
#'
#' @details Groups:
#' 2=c('CFYWMLIV', 'GPATSNHQEDRK'),
#' 3=c('CFYWMLIV', 'GPATS', 'NHQEDRK'),
#' 4=c('CFYW', 'MLIV', 'GPATS', 'NHQEDRK'),
#' 5=c('CFYW', 'MLIV', 'G', 'PATS', 'NHQEDRK'),
#' 6=c('CFYW', 'MLIV', 'G', 'P', 'ATS', 'NHQEDRK'),
#' 7=c('CFYW', 'MLIV', 'G', 'P', 'ATS', 'NHQED', 'RK'),
#' 8=c('CFYW', 'MLIV', 'G', 'P', 'ATS', 'NH', 'QED', 'RK'),
#' 9=c('CFYW', 'ML', 'IV', 'G', 'P', 'ATS', 'NH', 'QED', 'RK'),
#' 10=c('C', 'FYW', 'ML', 'IV', 'G', 'P', 'ATS', 'NH', 'QED', 'RK'),
#' 11=c('C', 'FYW', 'ML', 'IV', 'G', 'P', 'A', 'TS', 'NH', 'QED', 'RK'),
#' 12=c('C', 'FYW', 'ML', 'IV', 'G', 'P', 'A', 'TS', 'NH', 'QE', 'D', 'RK'),
#' 13=c('C', 'FYW', 'ML', 'IV', 'G', 'P', 'A', 'T', 'S', 'NH', 'QE', 'D', 'RK'),
#' 14=c('C', 'FYW', 'ML', 'IV', 'G', 'P', 'A', 'T', 'S', 'N', 'H', 'QE', 'D', 'RK'),
#' 15=c('C', 'FYW', 'ML', 'IV', 'G', 'P', 'A', 'T', 'S', 'N', 'H', 'QE', 'D', 'R', 'K'),
#' 16=c('C', 'FY', 'W', 'ML', 'IV', 'G', 'P', 'A', 'T', 'S', 'N', 'H', 'QE', 'D', 'R', 'K'),
#' 17=c('C', 'FY', 'W', 'ML', 'IV', 'G', 'P', 'A', 'T', 'S', 'N', 'H', 'Q', 'E', 'D', 'R', 'K'),
#' 18=c('C', 'FY', 'W', 'M', 'L', 'IV', 'G', 'P', 'A', 'T', 'S', 'N', 'H', 'Q', 'E', 'D', 'R', 'K'),
#' 19=c('C', 'F', 'Y', 'W', 'M', 'L', 'IV', 'G', 'P', 'A', 'T', 'S', 'N', 'H', 'Q', 'E', 'D', 'R', 'K'),
#' 20=c('C', 'F', 'Y', 'W', 'M', 'L', 'I', 'V', 'G', 'P', 'A', 'T', 'S', 'N', 'H', 'Q', 'E', 'D', 'R', 'K')
#'
#' @export
#'
#' @examples
#'
#' filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
#'
#' mat1<-PseKRAAC_T11(seqs=filePrs,type="gap",Grp=4,GapOrLambdaValue=3,k=2)
#'
#' mat2<-PseKRAAC_T11(seqs=filePrs,type="lambda",Grp=4,GapOrLambdaValue=3,k=2)


PseKRAAC_T11 <- function(seqs,type="gap",Grp=5,GapOrLambdaValue=2,k=4,label=c())
{


  numGrp=Grp

  if(length(seqs)==1&&file.exists(seqs)){
    seqs<-fa.read(seqs,alphabet="aa")
    seqs_Lab<-alphabetCheck(seqs,alphabet = "aa",label)

    seqs<-seqs_Lab[[1]]
    label<-seqs_Lab[[2]]
  }
  else if(is.vector(seqs)){
    seqs<-sapply(seqs,toupper)

    seqs_Lab<-alphabetCheck(seqs,alphabet = "aa",label)

    seqs<-seqs_Lab[[1]]
    label<-seqs_Lab[[2]]

  }
  else {
    stop("ERROR: Input sequence is not in the correct format. It should be a FASTA file or a string vector.")
  }


  if(type!="lambda"&&type!="gap"){
    stop("ERROR: type should be lambda or gap ")
  }

  lenSeqs<-sapply(seqs,nchar)
  numSeqs<-length(seqs)

  Grp<-paste0("Grp",Grp)

  AAGroup = list(
    Grp2=c('CFYWMLIV', 'GPATSNHQEDRK'),
    Grp3=c('CFYWMLIV', 'GPATS', 'NHQEDRK'),
    Grp4=c('CFYW', 'MLIV', 'GPATS', 'NHQEDRK'),
    Grp5=c('CFYW', 'MLIV', 'G', 'PATS', 'NHQEDRK'),
    Grp6=c('CFYW', 'MLIV', 'G', 'P', 'ATS', 'NHQEDRK'),
    Grp7=c('CFYW', 'MLIV', 'G', 'P', 'ATS', 'NHQED', 'RK'),
    Grp8=c('CFYW', 'MLIV', 'G', 'P', 'ATS', 'NH', 'QED', 'RK'),
    Grp9=c('CFYW', 'ML', 'IV', 'G', 'P', 'ATS', 'NH', 'QED', 'RK'),
    Grp10=c('C', 'FYW', 'ML', 'IV', 'G', 'P', 'ATS', 'NH', 'QED', 'RK'),
    Grp11=c('C', 'FYW', 'ML', 'IV', 'G', 'P', 'A', 'TS', 'NH', 'QED', 'RK'),
    Grp12=c('C', 'FYW', 'ML', 'IV', 'G', 'P', 'A', 'TS', 'NH', 'QE', 'D', 'RK'),
    Grp13=c('C', 'FYW', 'ML', 'IV', 'G', 'P', 'A', 'T', 'S', 'NH', 'QE', 'D', 'RK'),
    Grp14=c('C', 'FYW', 'ML', 'IV', 'G', 'P', 'A', 'T', 'S', 'N', 'H', 'QE', 'D', 'RK'),
    Grp15=c('C', 'FYW', 'ML', 'IV', 'G', 'P', 'A', 'T', 'S', 'N', 'H', 'QE', 'D', 'R', 'K'),
    Grp16=c('C', 'FY', 'W', 'ML', 'IV', 'G', 'P', 'A', 'T', 'S', 'N', 'H', 'QE', 'D', 'R', 'K'),
    Grp17=c('C', 'FY', 'W', 'ML', 'IV', 'G', 'P', 'A', 'T', 'S', 'N', 'H', 'Q', 'E', 'D', 'R', 'K'),
    Grp18=c('C', 'FY', 'W', 'M', 'L', 'IV', 'G', 'P', 'A', 'T', 'S', 'N', 'H', 'Q', 'E', 'D', 'R', 'K'),
    Grp19=c('C', 'F', 'Y', 'W', 'M', 'L', 'IV', 'G', 'P', 'A', 'T', 'S', 'N', 'H', 'Q', 'E', 'D', 'R', 'K'),
    Grp20=c('C', 'F', 'Y', 'W', 'M', 'L', 'I', 'V', 'G', 'P', 'A', 'T', 'S', 'N', 'H', 'Q', 'E', 'D', 'R', 'K')
  )

  if(is.null(AAGroup[[Grp]])){
    stop(paste0("Error: There is no ",Grp))
  }

  lenGrp<-length(AAGroup[[Grp]])
  aa<-vector()
  chars<-vector()
  VectGrp<-c("Grp10"='a',"Grp11"='b',"Grp12"='c',"Grp13"='d',"Grp14"='e',"Grp15"='f',"Grp16"='g',"Grp17"='h',"Grp18"='i',"Grp19"='j',"Grp20"='k')
  for(i in 1:lenGrp){
    str<-AAGroup[[Grp]][i]
    lenStr<-nchar(str)
    if(i<10){
      aa<-c(aa,rep(i,lenStr))
    } else {
      aa<-c(aa,rep(VectGrp[(i-9)],lenStr))
    }
    chars<-c(chars,unlist(strsplit(str,split = "")))
  }

  names(aa)<-chars

  grpChars<-lapply(seqs,function(s){
    Chars<-unlist(strsplit(s,split = ""))
    return(aa[Chars])
  })



  if(type=="lambda"){
    lambda_model(grpchars = grpChars,GapOrLambdaValue = GapOrLambdaValue,k = k,numGrp = numGrp)
  }
  else if(type=="gap"){
    grpStr<-lapply(grpChars, function(s) paste0(s,collapse = ""))
    gap_model(grpstr  = grpStr,GapOrLambdaValue = GapOrLambdaValue,k = k,numGrp = numGrp)
  }
  else{
    stop("type should be lambda or gap")
  }



}




lambda_model<-function(grpchars,GapOrLambdaValue,k,numGrp){

  featureMatrix<-matrix(0,ncol = (numGrp^k),nrow = length(grpchars))
  colnames(featureMatrix)<-nameKmer(k=k,type = "num",num=numGrp)
  lens<-lapply(grpchars, length)
  strIndx=1
  for(n in 1:length(grpchars)){
    len<-lens[n]
    endIdx<-strIndx+((k-1)*(GapOrLambdaValue+1))
    if(endIdx>len){
      stop(paste0("ERROR: Sequence ",n," doesn't have a suffitient length"))
    }
    tempVect<-seq(1,endIdx,(GapOrLambdaValue+1))

    while(endIdx<=len){
      name<-paste(grpchars[[n]][tempVect],collapse = '')
      featureMatrix[n,name]<-featureMatrix[n,name]+1
      tempVect<-tempVect+1
      endIdx<-endIdx+1
    }

  }
  return(featureMatrix)
}

gap_model<-function(grpstr,GapOrLambdaValue,k,numGrp){

  featureMatrix<-matrix(0,ncol = (numGrp^k),nrow = length(grpstr))
  colnames(featureMatrix)<-nameKmer(k=k,type = "num",num=numGrp)

  substrVector<-lapply(grpstr, function(x){
    numChars<-nchar(x)
    if(numChars<k){
      stop(paste0("ERROR: Length of all the sequences should be greater than k"))
    }
    strts<-(numChars-k+1)
    g<-(GapOrLambdaValue+k)

    substrVector<-substring(x, seq(1, strts, g), seq(k, numChars, g))
    return(substrVector)
  })
  tabSubVect<-lapply(substrVector,table)
  for(j in 1:length(grpstr)){
    featureMatrix[j,names(tabSubVect[[j]])]<-tabSubVect[[j]]

  }


  return(featureMatrix)

}
