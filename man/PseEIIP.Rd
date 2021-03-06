% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PseEIIP.R
\name{PseEIIP}
\alias{PseEIIP}
\title{Pseudo Electron-Ion Interaction Pseudopotentials of Trinucleotide (PseEIIP)}
\usage{
PseEIIP(seqs, label = c())
}
\arguments{
\item{seqs}{is a FASTA file containing nucleotide sequences. The sequences start
with '>'. Also, seqs could be a string vector. Each element of the vector is a nucleotide sequence.}

\item{label}{is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
each entry (i.e., sequence).}
}
\value{
This function returns a feature matrix which the number of rows is equal to the number of sequences and the
number of columns is 4^3=64.
}
\description{
This function calculates the pseudo electron-ion interaction for each sequence.
It creates a feature vector for each sequence.
The vector contains a value for each for each tri-nucleotide.
The value is computed by multiplying the aggregate value of electron-ion interaction of each nucleotide
}
\examples{

LNCSeqsADR<-system.file("extdata/",package="ftrCOOL")
LNC50Nuc<-as.vector(read.csv(paste0(LNCSeqsADR,"/LNC50Nuc.csv"))[,2])
mat<-PseEIIP(seqs = LNC50Nuc)
}
\references{
Chen, Zhen, et al. "iLearn: an integrated platform and meta-learner for feature engineering, machine-learning analysis and modeling of DNA, RNA and protein sequence data." Briefings in bioinformatics 21.3 (2020): 1047-1057.
}
