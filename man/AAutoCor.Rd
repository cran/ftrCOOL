% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AAutoCor.R
\name{AAutoCor}
\alias{AAutoCor}
\title{Amino Acid Autocorrelation-Autocovariance (AAutoCor)}
\usage{
AAutoCor(
  seqs,
  selectedAAidx = list(c("CIDH920105", "BHAR880101", "CHAM820101", "CHAM820102",
    "CHOC760101", "BIGC670101", "CHAM810101", "DAYM780201")),
  maxlag = 3,
  threshold = 1,
  type = c("Moran", "Geary", "NormalizeMBorto", "AC", "CC", "ACC"),
  label = c()
)
}
\arguments{
\item{seqs}{is a FASTA file with amino acid sequences. Each sequence starts
with a '>' character. Also, seqs could be a string vector. Each element of the vector is a peptide/protein sequence.}

\item{selectedAAidx}{Function takes as input the physicochemical properties. Users select the properties by their ids
or indices in the aaIndex2 file. This parameter could be a vector or a list of amino acid indices.
The default values of the vector are the 'CIDH920105','BHAR880101','CHAM820101','CHAM820102','CHOC760101','BIGC670101','CHAM810101','DAYM780201'
ids in the aaIndex2 file.}

\item{maxlag}{This parameter shows the maximum gap between two amino acids. The gaps change from 1 to maxlag (the maximum lag).}

\item{threshold}{is a number between (0 , 1]. In selectedAAidx, indices with a correlation higher than the threshold will be deleted.
The default value is 1.}

\item{type}{could be 'Moran', 'Greay', 'NormalizeMBorto', 'AC', 'CC', or 'ACC'. Also, it could be any combination of them.}

\item{label}{is an optional parameter. It is a vector whose length is equivalent to the number of sequences. It shows the class of
each entry (i.e., sequence).}
}
\value{
This function returns a feature matrix. The number of columns in the matrix changes depending on the chosen autocorrelation or autocovariance types and nlag parameter.
The output is a matrix. The number of rows shows the number of sequences.
}
\description{
It creates the feature matrix for each function in autocorelation
(i.e., Moran, Greay, NormalizeMBorto) or autocovariance (i.e., AC, CC,ACC). The user can
select any combination of the functions too. In this case, the final matrix will contain
features of each selected function.
}
\details{
For CC and AAC autocovriance functions, which consider the covariance of the
two physicochemical properties, we have provided users with the ability to categorize
their selected properties in a list.
The binary combination of each group will be taken into account.
Note: If all the features are in a group or selectedAAidx parameter is a vector,
the binary combination will be calculated for all the physicochemical properties.
}
\examples{

filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
mat1<-AAutoCor(seqs=filePrs,maxlag=20,threshold=0.9,
type=c("Moran","Geary","NormalizeMBorto","AC"))

mat2<-AAutoCor(seqs=filePrs,maxlag=20,threshold=0.9,selectedAAidx=
list(c('CIDH920105','BHAR880101','CHAM820101','CHAM820102'),c('CHOC760101','BIGC670101')
,c('CHAM810101','DAYM780201')),type=c("AC","CC","ACC"))

}
