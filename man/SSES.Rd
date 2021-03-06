% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SSES.R
\name{SSES}
\alias{SSES}
\title{Secondary Structure Elements Simple (SSES)}
\usage{
SSES(dirPath, outFormat = "mat", outputFileDist = "")
}
\arguments{
\item{dirPath}{Path of the directory which contains all output files of PSIPRED. Each file belongs to a sequence.}

\item{outFormat}{It can take two values: 'mat' (which stands for matrix) and 'txt'. The default value is 'mat'.}

\item{outputFileDist}{It shows the path and name of the 'txt' output file.}
}
\value{
The output depends on the outFormat which can be either 'mat' or 'txt'. If outFormat is 'mat', the function returns a feature
matrix for sequences with the same lengths such that the number of columns is equal to the length of the sequences
and the number of rows is equal to the number of sequences.
If the outFormat is 'txt', the output is written to a tab-delimited file.
}
\description{
This function works based on the output of PSIPRED which predicts the secondary structure of the amino acids in a sequence.
The output of the PSIPRED is a tab-delimited file which contains the secondary structure in the third column.
The function represent amino acids in the helix structure by 'H', amino acids in the extended structure by 'E', and amino acids in the coil structure by 'C'.
}
\note{
This function is provided for the sequences with the same lengths. However,
the users can use 'txt' option in the outFormat parameter for sequences with different lengths.
Warning: If the outFormat is set to 'mat' for sequences with different lengths, it returns an error.
Also, when the output format is 'txt', the label information is not displayed in the text file.
It is noteworthy that, 'txt' format is not usable for machine learning purposes.
}
\examples{

dir = tempdir()
ad<-paste0(dir,"/simpleSSE.txt")

Predss2dir<-system.file("testForder",package="ftrCOOL")
Predss2dir<-paste0(Predss2dir,"/ss2Dir/")
mat<-SSES(Predss2dir,outFormat="txt",outputFileDist=ad)

unlink("dir", recursive = TRUE)
}
