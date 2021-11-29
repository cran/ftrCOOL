#' Multivariate Mutual Information_DNA (MMI_DNA)
#'
#' MMI computes mutual information based on 2-mers T2 = { AA, AC, AG, AT, CC, CG, CT, GG, GT, TT}
#' and 3-mers T3 = {AAA, AAC, AAG, AAT, ACC, ACG, ACT, AGG, AGT, ATT, CCC, CCG, CCT, CGG, CGT, CTT, GGG, GGT, GTT and TTT}
#' for more information please check the reference part.
#'
#' @references Zhen Chen, Pei Zhao, Chen Li, Fuyi Li, Dongxu Xiang, Yong-Zi Chen, Tatsuya Akutsu, Roger J Daly, Geoffrey I Webb, Quanzhi Zhao, Lukasz Kurgan, Jiangning Song. iLearnPlus: a comprehensive and automated machine-learning platform for nucleic acid and protein sequence analysis, prediction and visualization, Nucleic Acids Research (2021).
#'
#'
#'
#' @param seqs is a FASTA file containing nucleotide sequences. The sequences start
#' with '>'. Also, seqs could be a string vector. Each element of the vector is a nucleotide sequence.
#'
#'
#' @param label is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
#' each entry (i.e., sequence).
#'
#'
#'
#' @return It is a feature matrix. The number of columns is 30 and the number of rows is equal to the number of sequences.
#'
#' @export
#'
#'
#' @examples
#'
#' fileLNC<-system.file("extdata/Athaliana_LNCRNA.fa",package="ftrCOOL")
#' mat<-MMI_DNA(seqs=fileLNC)




MMI_DNA<- function(seqs,label=c())
{

  if(length(seqs)==1&&file.exists(seqs)){
    seqs<-fa.read(seqs,alphabet="dna")
    seqs_Lab<-alphabetCheck(seqs,alphabet = "dna",label)

    seqs<-seqs_Lab[[1]]
    label<-seqs_Lab[[2]]
  }
  else if(is.vector(seqs)){
    seqs<-sapply(seqs,toupper)

    seqs_Lab<-alphabetCheck(seqs,alphabet = "dna",label)

    seqs<-seqs_Lab[[1]]
    label<-seqs_Lab[[2]]

  }
  else {
    stop("ERROR: Input sequence is not in the correct format. It should be a FASTA file or a string vector.")
  }
  lenSeqs<-sapply(seqs,nchar)
  numSeqs=length(seqs)

  featureMatrix<-matrix(0, nrow = numSeqs, ncol = (10+20))

   diNames<-nameKmer(k=2,type = "dna")
   triNames<-nameKmer(k=3,type = "dna")



  DimerList<-list("AA"=c("AA"),"AC"=c("AC","CA"), "AG"=c("AG","GA"), "AT"=c("AT","TA"), "CC"=c("CC"), "CG"=c("CG","GC"), "CT"=c("CT","TC"), "GG"=c("GG"), "GT"=c("GT","TG"), "TT"=c("TT"))


  TrimerList<-list("AAA"=c("AAA"), "AAC"=c("AAC","ACA","CAA"), "AAG"=c("AAG","AGA","GAA"), "AAT"=c("AAT","ATA","TAA"), "ACC"=c("ACC","CAC","CCA"),"ACG"=c("ACG","AGC","CAG","CGA","GAC","GCA"), "ACT"=c("ACT","ATC","CAT","CTA","TAC","TCA"), "AGG"=c("AGG","GAG","GGA"), "AGT"=c("AGT","ATG","GAT","GTA","TAG","TGA"), "ATT"=c("ATT","TAT","TTA"), "CCC"=c("CCC"), "CCG"=c("CCG","CGC","GCC"), "CCT"=c("CCT","CTC","TCC"), "CGG"=c("CGG","GCG","GGC"), "CGT"=c("CGT","CTG","GCT","GTC","TCG","TGC"), "CTT"=c("CTT","TCT","TTC"), "GGG"=c("GGG"), "GGT"=c("GGT","GTG","TGG"), "GTT"=c("GTT","TGT","TTG"), "TTT"=c("TTT"))


  for(n in 1:numSeqs){
    seq<-seqs[n]
    lenSeq<-lenSeqs[n]
    chars<-unlist(strsplit(seq,split = ""))
    tabChars<-table(chars)/lenSeq

    temp1<-chars[1:(lenSeq-1)]
    temp2<-chars[2:lenSeq]
    Dimer<-paste0(temp1,temp2)
    tabDimer<-table(Dimer)/(lenSeq-1)
    wholeVectDi<-vector(mode = "numeric",16)
    names(wholeVectDi)<-diNames
    wholeVectDi[names(tabDimer)]<-tabDimer

    newVectDi<-vector(mode = "numeric",length = 10)
    names(newVectDi)<-names(DimerList)

    for(i in 1:10){
      x=wholeVectDi[DimerList[[i]]]
      newVectDi[i]<-sum(x)
    }


    Idi<-vector()
    for(i in 1:10){
      Dinuc<-unlist(strsplit(names(newVectDi)[i],split = ""))
      if(newVectDi[i]==0){
        var<-0
        Idi<-c(Idi,var)
        next()
      }
      temp<-log((newVectDi[i]/(tabChars[Dinuc[1]]*tabChars[Dinuc[2]])))
      var<-newVectDi[i]*temp
      Idi<-c(Idi,var)
    }

    names(Idi)<-names(newVectDi)

    temp1<-chars[1:(lenSeq-2)]
    temp2<-chars[2:(lenSeq-1)]
    temp3<-chars[3:lenSeq]
    Trimer<-paste0(temp1,temp2,temp3)
    tabTrimer<-table(Trimer)/(lenSeq-2)
    wholeVectTri<-vector(mode = "numeric",64)

    names(wholeVectTri)<-triNames
    wholeVectTri[names(tabTrimer)]<-tabTrimer

    newVectTri<-vector(mode = "numeric",length = 20)
    names(newVectTri)<-names(TrimerList)

    for(i in 1:20){
      x=wholeVectTri[TrimerList[[i]]]
      newVectTri[i]<-sum(x)
    }



    Itri<-vector()


    for(i in 1:20){

      trinuc<-unlist(strsplit(names(TrimerList)[i],split = ""))


      N1N2N3<-paste0(trinuc[1],trinuc[2],trinuc[3])
      fN1N2N3<-newVectTri[N1N2N3]

      N1N2<-paste0(trinuc[1],trinuc[2])
      fN1N2<-newVectDi[N1N2]

      N1N3<-paste0(trinuc[1],trinuc[3])
      fN1N3<-newVectDi[N1N3]

      N2N3<-paste0(trinuc[2],trinuc[3])
      fN2N3<-newVectDi[N2N3]

      fN1<-tabChars[trinuc[1]]

      fN2<-tabChars[trinuc[2]]

      fN3<-tabChars[trinuc[3]]

      var1<-fN1N2*log(fN1N2/(fN1*fN2))

      if(fN1N2==0)
        var1<-0
      var2<-(fN1N3/fN3)*log(fN1N3/fN3)
      if(fN1N3==0)
        var2<-0
      var3<-(fN1N2N3/fN2N3)*log(fN1N2N3/fN2N3)
      if(fN1N2N3==0)
        var3<-0


      Itri<-c(Itri,var1+var2-var3)
    }


    MMI<-c(Idi,Itri)
    featureMatrix[n,]<-MMI

  }
  Nam<-c(names(DimerList),names(TrimerList))
  temp<-paste0("MMI_",Nam)
  colnames(featureMatrix)<-temp

  if(length(label)==numSeqs){
    featureMatrix<-as.data.frame(featureMatrix)
    featureMatrix<-cbind(featureMatrix,label)
  }
  row.names(featureMatrix)<-names(seqs)
  return(featureMatrix)

}

