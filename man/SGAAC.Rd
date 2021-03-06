% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SGAAC.R
\name{SGAAC}
\alias{SGAAC}
\title{Splitted Group Amino Acid Composition (SGAAC)}
\usage{
SGAAC(
  seqs,
  k = 1,
  numNterm = 25,
  numCterm = 25,
  Grp = "locFus",
  normalized = TRUE,
  label = c()
)
}
\arguments{
\item{seqs}{is a FASTA file with amino acid sequences. Each sequence starts
with a '>' character. Also, seqs could be a string vector. Each element of the vector is a peptide/protein sequence.}

\item{k}{shows which type of amino acid composition applies to the parts.
For example, the amino acid composition is applied when k=1 and when k=2, the dipeptide Composition is applied.}

\item{numNterm}{shows how many amino acids should be considered for N-terminal.}

\item{numCterm}{shows how many amino acids should be considered for C-terminal.}

\item{Grp}{is a list of vectors containig amino acids. Each vector represents a category. Users can define a customized amino acid grouping, provided that the sum of all amino acids is 20 and there is no repeated amino acid in the groups.
Also, users can choose 'cTriad'(conjointTriad), 'locFus', or 'aromatic'. Each option provides specific information about the type of an amino acid grouping.}

\item{normalized}{is a logical parameter. When it is FALSE, the return value of the function does not change. Otherwise, the return value is normalized using the length of the sequence.}

\item{label}{is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
each entry (i.e., sequence).}
}
\value{
It returns a feature matrix. The number of rows is equal to the number of sequences.
The number of columns is 3*((number of groups)^k).
}
\description{
In this function, amino acids are first grouped into a user-defined category.
Later, the splitted amino Acid composition is computed.
Please note that this function differs from \link{SAAC} which works on individual amino acids.
}
\examples{
filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
mat<-SGAAC(seqs=filePrs,k=1,numNterm=15,numCterm=15,Grp="aromatic")

}
