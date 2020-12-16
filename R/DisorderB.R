#' disorder Binary (DisorderB)
#'
#' This function extracts the ordered and disordered amino acids in protein or peptide sequences. The input to the function is provided by VSL2 software. Also, the function
#' converts order amino acids to '10' and disorder amino acids to '01'.
#'
#' @param binaryType It can take any of the following values: ('strBin','logicBin','numBin').
#' 'strBin' (String binary): each dinucleotide is represented by a string containing 2 characters(0-1). order = "10"   disorder="01".
#' 'logicBin' (logical value): Each amino acid is represented by a vector containing 2 logical entries.  order = c(TRUE,FALSE)   disorder=c(FALSE,TRUE).
#' 'numBin' (numeric bin): Each amino acid is represented by a numeric (i.e., integer) vector containing 2 numeric entries.   order = c(1,0)   disorder=c(0,1).
#'
#' @note This function is provided for sequences with the same lengths.
#' Users can use 'txt' option in outFormat for sequences with different lengths.
#' Warning: If outFormat is set to 'mat' for sequences with different lengths, it returns an error.
#' Also, when output format is 'txt', label information is not shown in the text file.
#' It is noteworthy that 'txt' format is not usable for machine learning purposes if sequences have different sizes. Otherwise 'txt' format
#' is also usable for machine learning purposes.
#'
#' @param dirPath Path of the directory which contains all output files of VSL2. Each file belongs to a sequence.
#'
#' @param outFormat It can take two values: 'mat' (which stands for matrix) and 'txt'. The default value is 'mat'.
#'
#' @param outputFileDist It shows the path and name of the 'txt' output file.
#'
#'
#' @return The output is different depending on the outFormat parameter ('mat' or 'txt').
#' If outFormat is set to 'mat', it returns a feature matrix for sequences with the same lengths.
#' The number of rows is equal to the number of sequences and if binaryType is 'strBin', the number of columns is the length of the sequences.
#' Otherwise, it is equal to (length of the sequences)*2.
#' If outFormat is 'txt', all binary values will be written to a tab-delimited file. Each line in the file shows the binary format of a sequence.
#'
#'
#' @export
#'
#' @examples
#'
#' dir = tempdir()
#' PredDisdir<-system.file("testForder",package="ftrCOOL")
#'
#' PredDisdir<-paste0(PredDisdir,"/Disdir/")
#' ad1<-paste0(dir,"/disorderB.txt")
#'
#' mat<-DisorderB(PredDisdir,binaryType="strBin",outFormat="txt",outputFileDist=ad1)
#'
#'
DisorderB<-function(dirPath,binaryType="strBin",outFormat="mat",outputFileDist=""){


  disorderVectorSimple<-readDisDir(dirPath)
  lens<-lapply(disorderVectorSimple, length)

  #lens<-unlist(lens)

  if(outFormat=="mat"){

      if(length(unique(lens))>1){
        stop("ERROR: All sequences should have the same length in 'mat' mode. For sequences with different lengths, please use 'txt' for outFormat parameter")
      }

      if(binaryType=="strBin"){
        simpleMatrix=matrix("",ncol = lens[[1]],nrow = length(disorderVectorSimple))
        colnames(simpleMatrix)<-paste0("pos",1:lens[[1]])
        for(i in 1:nrow(simpleMatrix)){
          vect<-disorderVectorSimple[[i]]
          vect[vect=="D"]="01"
          vect[vect=="O"]="10"
          simpleMatrix[i,]=vect
        }
        return(simpleMatrix)
      }
      else if(binaryType=="logicalBin"){
        nucs=list("D"=c(FALSE,TRUE),"O"=c(TRUE,FALSE))

        featureList<-lapply(disorderVectorSimple,function(x) {
          cods<-nucs[x]
          cods<-unlist(cods)
          return(cods)
        })
        featureMatrix<-matrix(unlist(featureList),byrow = TRUE,ncol = (2*lens[[1]]),nrow = length(disorderVectorSimple))
        return(featureMatrix)

      } else if(binaryType=="numBin"){
        nucs=list("D"=c(0,1),"O"=c(1,0))

        featureList<-lapply(disorderVectorSimple,function(x) {
          cods<-nucs[x]
          cods<-unlist(cods)
          return(cods)
        })
        featureMatrix<-matrix(unlist(featureList),byrow = TRUE,ncol = (2*lens[[1]]),nrow = length(disorderVectorSimple))
        return(featureMatrix)
      }

    } else {
      for(i in 1:length(disorderVectorSimple)){
        vect<-disorderVectorSimple[[i]]
        vect[vect=="D"]="01"
        vect[vect=="O"]="10"
        tempName<-paste0("sequence",i)
        vect<-c(tempName,vect)
        temp<-paste(vect,collapse = "\t")
        write(temp,outputFileDist,append = TRUE)
      }
    }


}

