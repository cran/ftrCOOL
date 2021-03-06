% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LocalPoSpKNUCF_RNA.R
\name{LocalPoSpKNUCF_RNA}
\alias{LocalPoSpKNUCF_RNA}
\title{Local Position Specific k riboNucleotide Frequency (LocalPoSpKNUCF_RNA)}
\usage{
LocalPoSpKNUCF_RNA(
  seqs,
  k = 2,
  label = c(),
  outFormat = "mat",
  outputFileDist = ""
)
}
\arguments{
\item{seqs}{is a FASTA file containing ribonucleotide sequences. The sequences start
with '>'. Also, seqs could be a string vector. Each element of the vector is a ribonucleotide sequence.}

\item{k}{is a numeric value which holds the value of k in the k-mers.}

\item{label}{is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
each entry (i.e., sequence).}

\item{outFormat}{(output format) can take two values: 'mat'(matrix) and 'txt'. The default value is 'mat'.}

\item{outputFileDist}{shows the path and name of the 'txt' output file.}
}
\value{
The output depends on the outFormat parameter which can be either 'mat' or 'txt'. If outFormat is 'mat', the function returns a feature
matrix for sequences with the same length such that the number of columns is (sequence length-k+1)
and the number of rows is equal to the number of sequences.
If the outFormat is 'txt', the output is written to a tab-delimited file.
}
\description{
For each sequence, this function creates a feature vector denoted as
(f1,f2, f3, …, fN), where fi = freq(i'th k-mer of the sequence) / i.
It should be applied to sequences with the same length.
}
\note{
This function is provided for sequences with the same lengths.
Users can use 'txt' option in outFormat for sequences with different lengths.
Warning: If outFormat is set to 'mat' for sequences with different lengths, it returns an error.
Also, when output format is 'txt', label information is not shown in the text file.
It is noteworthy that 'txt' format is not usable for machine learning purposes if sequences have different sizes. Otherwise 'txt' format
is also usable for machine learning purposes.
}
\examples{

dir = tempdir()
fileLNC<-system.file("extdata/Carica_papaya101RNA.txt",package="ftrCOOL")
mat<-LocalPoSpKNUCF_RNA(seqs = fileLNC, k=2,outFormat="mat")

ad<-paste0(dir,"/LocalPoSpKnucF.txt")
fileLNC<-system.file("extdata/Carica_papaya101RNA.txt",package="ftrCOOL")
LocalPoSpKNUCF_RNA(seqs = fileLNC,k=1,outFormat="txt"
,outputFileDist=ad)

unlink("dir", recursive = TRUE)
}
