% --- Source file: subject.list.Rd ---
\name{subject.list}
\alias{subject.list}
\docType{data}
\title{
 List to describe the file of subject ids 
}
\description{
 The list to describe the file of subject ids for \code{\link{geno.list}}
}

\format{
        The format is:
         List of 4
         \describe{
          \item{file}{ Text file containing the subject ids.
             The file can be a single column of ids, or a delimited file
                   of several columns with the ids as one of the columns. No default.}
          \item{id.var}{ Column number(s) or variable name(s) containing the subject ids.
                         Use \code{id.var}=-1 if the file is a single column of ids.}
          \item{delimiter}{ The delimiter in \code{file}.
                         The default is "".}
         \item{header}{0 or 1 if the file contains a header of variable names. 
                       The default is 0.}
         }
}


\details{
  The order of the ids in this file must match the order of the genotypes in the genotype file.
  When using the \code{\link{scan.score}} function, this list can often just be set to the name of the file
  containing the ids.
}

\keyword{misc}


