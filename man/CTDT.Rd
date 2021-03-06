% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CTDT.R
\name{CTDT}
\alias{CTDT}
\title{CTD Transition (CTDT)}
\usage{
CTDT(seqs, normalized = FALSE, label = c())
}
\arguments{
\item{seqs}{is a FASTA file with amino acid sequences. Each sequence starts
with a '>' character. Also, seqs could be a string vector. Each element of the vector is a peptide/protein sequence.}

\item{normalized}{is a logical parameter. When it is FALSE, the return value of the function does not change. Otherwise, the return value is normalized using the length of the sequence.}

\item{label}{is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
each entry (i.e., sequence).}
}
\value{
This function returns a feature matrix.
The number of rows is equal to the number of sequences and the number of columns is 3*7,
where three is the number of transition types (i.e., 1-2, 1-3, and 2-3) and thirteen is the number of properties.
}
\description{
This function computes the transition part of \link{CTD}. Thirteen properties are defined in this function. Each
property categorizes the amino acids of a sequence into three groups. For each property, the grouped amino acid transition (i.e., transitions 1-2, 1-3, and 2-3)
is calculated. For more information, please check the references.
}
\examples{

filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
CTD_T<-CTDT(seqs=filePrs,normalized=FALSE)

}
\references{
Dubchak, Inna, et al. "Prediction of protein folding class using global description of amino acid sequence." Proceedings of the National Academy of Sciences 92.19 (1995): 8700-8704.
}
