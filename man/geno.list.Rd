% --- Source file: geno.list.Rd ---
\name{geno.list}
\alias{geno.list}
\docType{data}
\title{
 List to describe the genotype data
}
\description{
 The list to describe the genotype data for \code{\link{scan.score}}.
}
\format{
      The format is: List of 4
       \describe{
           \item{file}{File to use. Currently, only output files from IMPUTE2 
               containing genotype probabilities can be processed. 
                 No default. }
           \item{subject.list}{ List to describe the subject ids stored in a file.
                             The order of the subject ids is \code{subject.list$file}
                             must match the order in the genotype data. See \code{\link{subject.list}}.
                             No default.}
           \item{start.vec}{Starting row of \code{file} to begin processing SNPs. The default is 1. }
           \item{stop.vec}{Last row of \code{file} to finish processing SNPs. 
                           Use any value < 1 or \code{Inf} so that all SNPs from rows \code{start.vec} 
                           to the end of the file will be analyzed.
                           The default is -1. }
       }
}


\details{ In this list, \code{file} and \code{subject.list} must be specified, 
}
%\examples{}
\keyword{misc}

