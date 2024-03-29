% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MMI_RNA.R
\name{MMI_RNA}
\alias{MMI_RNA}
\title{Multivariate Mutual Information_RNA (MMI_RNA)}
\usage{
MMI_RNA(seqs, label = c())
}
\arguments{
\item{seqs}{is a FASTA file containing ribonucleotide sequences. The sequences start
with '>'. Also, seqs could be a string vector. Each element of the vector is a ribonucleotide sequence.}

\item{label}{is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
each entry (i.e., sequence).}
}
\value{
It is a feature matrix. The number of columns is 30 and the number of rows is equal to the number of sequences.
}
\description{
MMI computes mutual information based on 2-mers T2 = { AA, AC, AG, AU, CC, CG, CU, GG, GU, U}
and 3-mers T3 = {AAA, AAC, AAG, AAU, ACC, ACG, ACU, AGG, AGU, AUU, CCC, CCG, CCU, CGG, CGU, CUU, GGG, GGU, GUU and UUU}
for more information please check the reference part.
}
\examples{

fileLNC<-system.file("extdata/Carica_papaya101RNA.txt",package="ftrCOOL")
mat<-MMI_RNA(seqs=fileLNC)
}
\references{
Zhen Chen, Pei Zhao, Chen Li, Fuyi Li, Dongxu Xiang, Yong-Zi Chen, Tatsuya Akutsu, Roger J Daly, Geoffrey I Webb, Quanzhi Zhao, Lukasz Kurgan, Jiangning Song. iLearnPlus: a comprehensive and automated machine-learning platform for ribonucleic acid and protein sequence analysis, prediction and visualization, Nucleic Acids Research (2021).
}
