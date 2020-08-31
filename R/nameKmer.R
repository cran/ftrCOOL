#' nameKmer
#'
#' This function creates all possible k-combinations of the given alphabets.
#'
#'
#' @param k is a numeric value.
#'
#' @param type can be one of "aa", "rna", "dna", or "num".
#'
#' @param num When type is set to "num", it shows the numeric alphabet( 1,..,,num).
#'
#' @return a string vector of length (20^k for 'aa' type), (4^k for 'dna' type), (4^k for 'rna' type),
#' and (num^k for 'num' type).
#'
#'
#' @export
#'
#' @examples
#'
#' all_kmersAA<-nameKmer(k=2,type="aa")
#'
#' all_kmersDNA<-nameKmer(k=3,type="dna")
#'
#' all_kmersNUM<-nameKmer(k=3,type="num",num=2)

nameKmer<-function(k=3,type="aa",num=0){

  if(type=="aa"){
    vect<-c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
    len<-length(vect)
  } else if(type=="dna"){
    vect<-c("A","C","G","T")
    len<-length(vect)
  } else if(type=="rna"){
    vect<-c("A","C","G","T")
    len<-length(vect)
  } else if(type=="num") {
    if(num<10){
      vect<-1:num
      vect<-as.character(vect)
      len<-num
    } else {
      if(num>20){
        stop("ERROR: num should be less than 20")
      } else {
        VectGrp<-c("Grp10"='a',"Grp11"='b',"Grp12"='c',"Grp13"='d',"Grp14"='e',"Grp15"='f',"Grp16"='g',"Grp17"='h',"Grp18"='i',"Grp19"='j',"Grp20"='k')
        VectGrp<-VectGrp[1:(num-9)]
        vect<-c(1:9,VectGrp)
        len<-length(vect)
      }
    }


  }
  else {
    stop("ERROR: type should be 'aa'(amino acid) or 'rna' or 'dna'")
  }

  totalVect<-rep(vect,k)
  kmer<-sort(apply(expand.grid(split(totalVect, rep(1:k,each = len))), 1, paste, collapse = ""))
  return(kmer)

}
