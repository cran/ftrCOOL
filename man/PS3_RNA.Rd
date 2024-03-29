% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PS3_RNA.R
\name{PS3_RNA}
\alias{PS3_RNA}
\title{Position-specific of three ribonucleotide_RNA (PS3_RNA)}
\usage{
PS3_RNA(
  seqs,
  binaryType = "numBin",
  label = c(),
  outFormat = "mat",
  outputFileDist = ""
)
}
\arguments{
\item{seqs}{is a FASTA file containing ribonucleotide sequences. The sequences start
with '>'. Also, seqs could be a string vector. Each element of the vector is a ribonucleotide sequence.}

\item{binaryType}{It can take any of the following values: ('strBin','logicBin','numBin').
'strBin'(String binary): each di-ribonucleotide is represented by a string containing 64 characters (63 times '0' and one '1'). For example, 'AAA' = "1000000000000000...0", ....
'logicBin'(logical value): Each amino acid is represented by a vector containing 64 logical entries (63 times F and one T). For example, 'AA' = c(T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,...,F), ...
'numBin' (numeric bin): Each amino acid is represented by a numeric (i.e., integer) vector containing 64 numerals (63 times '0' and one '1'). For example, 'AA' = c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,...,0)}

\item{label}{is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
each entry (i.e., sequence).}

\item{outFormat}{(output format) can take two values: 'mat'(matrix) and 'txt'. The default value is 'mat'.}

\item{outputFileDist}{shows the path and name of the 'txt' output file.}
}
\value{
The output is different depending on the outFormat parameter ('mat' or 'txt').
If outFormat is set to 'mat', it returns a feature matrix for sequences with the same lengths.
The number of rows is equal to the number of sequences and if binaryType is 'strBin', the number of columns is the length of the sequences.
Otherwise, it is equal to (length of the sequences-2)*64.
If outFormat is 'txt', all binary values will be written to a the output is written to a tab-delimited file. Each line in the file shows the binary format of a sequence.
}
\description{
This function transforms each tri-ribonucleotide of the sequence to a binary format.
The type of the binary format is determined by the binaryType parameter.
For details about each format, please refer to the description of the binaryType parameter.
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

fileLNC<-system.file("extdata/Carica_papaya101RNA.txt",package="ftrCOOL")
mat<-PS3_RNA(seqs = fileLNC, binaryType="numBin",outFormat="mat")




}
\references{
Zhen Chen, Pei Zhao, Chen Li, Fuyi Li, Dongxu Xiang, Yong-Zi Chen, Tatsuya Akutsu, Roger J Daly, Geoffrey I Webb, Quanzhi Zhao, Lukasz Kurgan, Jiangning Song, iLearnPlus: a comprehensive and automated machine-learning platform for nucleic acid and protein sequence analysis, prediction and visualization, Nucleic Acids Research, (2021).
}
