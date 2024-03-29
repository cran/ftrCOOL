% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DistancePair.R
\name{DistancePair}
\alias{DistancePair}
\title{PseAAC of distance-pairs and reduced alphabet (DistancePair)}
\usage{
DistancePair(seqs, rng = 3, normalized = TRUE, Grp = "cp14", label = c())
}
\arguments{
\item{seqs}{is a FASTA file with amino acid sequences. Each sequence starts
with a '>' character. Also, seqs could be a string vector. Each element of the vector is a peptide/protein sequence.}

\item{rng}{This parameter is a number. It shows maximum number of spaces between amino acid pairs.
For each k in the rng vector, a new vector (whose size is (number of categorizes)^2) is created which contains the frequency of pairs with k gaps.}

\item{normalized}{is a logical parameter. When it is FALSE, the return value of the function does not change. Otherwise, the return value is normalized using the length of the sequence.}

\item{Grp}{for this parameter users can choose between these items: 'cp13', 'cp14', 'cp19', or 'cp20'.}

\item{label}{is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
each entry (i.e., sequence).}
}
\value{
This function returns a feature matrix. Row length is equal to the number of sequences and
the number of columns is (number of categorizes)+((number of categorizes)^2)*(rng+1).
}
\description{
In this function, first amino acids are grouped into a category which is one of 'cp13', 'cp14', 'cp19', 'cp20'. Users choose one of these terms to categorize amino acids.
Then DistancePair function computes frequencies of all grouped residues and also all grouped-paired residues with [0,rng] distance. 'rng'
is a parameter which already was set by the user.
}
\examples{

filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
mat1<-DistancePair(seqs=filePrs,rng=2,Grp="cp14")
}
