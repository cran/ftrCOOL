% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Mismatch_RNA.R
\name{Mismatch_RNA}
\alias{Mismatch_RNA}
\title{Mismatch_RNA (Mismatch_RNA)}
\usage{
Mismatch_RNA(seqs, k = 3, m = 2, label = c())
}
\arguments{
\item{seqs}{is a FASTA file containing ribonucleotide sequences. The sequences start
with '>'. Also, seqs could be a string vector. Each element of the vector is a ribonucleotide sequence.}

\item{k}{This parameter can be a number which shows kmer.}

\item{m}{This parametr shows muximum number of mismatches.}

\item{label}{is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
each entry (i.e., sequence).}
}
\value{
This function returns a feature matrix. The number of rows is equal to the number of sequences and
the number of columns depends on the rng vector. For each value k in the vector, (4)^k columns are created in the matrix.
}
\description{
This function also calculates the frequencies of all k-mers in the sequence but alows maximum m mismatch.
m<k.
}
\examples{

fileLNC<-system.file("extdata/Carica_papaya101RNA.txt",package="ftrCOOL")
mat<-Mismatch_RNA(seqs=fileLNC)


}
\references{
Liu, B., Gao, X. and Zhang, H. BioSeq-Analysis2.0: an updated platform for analyzing DNA, RNA and protein sequences at sequence level and residue level based on machine learning approaches. Nucleic Acids Res (2019).
}
