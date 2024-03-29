% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ENUComposition_DNA.R
\name{ENUComposition_DNA}
\alias{ENUComposition_DNA}
\title{Enhanced Nucleotide Composition (ENUComposition_DNA)}
\usage{
ENUComposition_DNA(
  seqs,
  winSize = 50,
  overLap = TRUE,
  label = c(),
  outFormat = "mat",
  outputFileDist = ""
)
}
\arguments{
\item{seqs}{is a FASTA file containing nucleotide sequences. The sequences start
with '>'. Also, seqs could be a string vector. Each element of the vector is a nucleotide sequence.}

\item{winSize}{is a number which shows the size of the window.}

\item{overLap}{This parameter shows how the window
moves on the sequence. If the overlap is set to TRUE, the next window would have distance 1 with
the previous window. Otherwise, the next window will start from the next nucleotide after the previous window.
There is no overlap between the next and previous windows.}

\item{label}{is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
each entry (i.e., sequence).}

\item{outFormat}{(output format) can take two values: 'mat'(matrix) and 'txt'. The default value is 'mat'.}

\item{outputFileDist}{shows the path and name of the 'txt' output file.}
}
\value{
The output depends on the outFormat parameter which can be either 'mat' or 'txt'. If outFormat is 'mat', the function returns a feature
matrix for sequences with the same length such that the number of columns is (4 * number of partitions displayed by the window)
and the number of rows is equal to the number of sequences.
If the outFormat is 'txt', the output is written to a tab-delimited file.
}
\description{
This function slides a window over the input sequence(s).
Also, it computes the composition of nucleotides that appears within the limits of the window.
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
LNCSeqsADR<-system.file("extdata/",package="ftrCOOL")
LNC50Nuc<-as.vector(read.csv(paste0(LNCSeqsADR,"/LNC50Nuc.csv"))[,2])
mat<-ENUComposition_DNA(seqs = LNC50Nuc, winSize=20,outFormat="mat")

ad<-paste0(dir,"/ENUCcompos.txt")
fileLNC<-system.file("extdata/Athaliana_LNCRNA.fa",package="ftrCOOL")
ENUComposition_DNA(seqs = fileLNC,outFormat="txt",winSize=20
,outputFileDist=ad,overLap=FALSE)


unlink("dir", recursive = TRUE)
}
