% --- Source file: pheno.list.Rd ---
\name{pheno.list}
\alias{pheno.list}
\docType{data}
\title{
 List to describe the covariate and outcome data 
}
\description{
 The list to describe the covariate and outcome data for \code{\link{scan.score}}.
}

\format{
        The format is:
         List of 7
         \describe{
          \item{file}{ Covariate data file. This file must have variable names that
                        include an id variable, case-control variable, and
                        subtype variables (see \code{id.var}, \code{cc.var} and \code{subtype.vars}).
                         No default.}
          \item{id.var}{ Name of the id variable(s). 
                         No default.}
          \item{cc.var}{ Name of the binary case-control variable. This variable
                        must be coded as 0 and 1, where 0 denotes a control.
                        No default. }
          \item{subtype.vars}{Character vector of the subtype variables. This vector
                              must have length greater than one. No Default.}
          \item{covars.obj}{NULL, a character vector of variable names, or a formula
                           for the adjusted covariates.
                        The default is NULL.}
          \item{delimiter}{ The delimiter in \code{file}.
                         The default is "".}
          \item{subsetData}{ List of sublists to subset the phenotype data for analyses. 
                             Each sublist should contain the names "var", "operator" and "value" corresponding
                             to a variable name, operator and values of the variable.
                             Multiple sublists are logically connected by the AND operator. For example, \cr
                          subsetData=list(list(var="GENDER", operator="==", value="MALE")) \cr
                          will only include subjects with the string "MALE" for the GENDER variable. \cr
                          subsetData=list(list(var="AGE", operator=">", value=50), \cr
                                          list(var="STUDY", operator="\%in\%", value=c("A", "B", "C"))) \cr
                          will include subjects with AGE > 50 AND in STUDY A, B or C.  
                          The default is NULL.}
         }
}


\details{ In this list, \code{file}, \code{id.var}, \code{cc.var}, and
   \code{subtype.vars} must be specified.
   The variable \code{id.var} is the link between the covariate data and the genotype data.
   For each subject id, there must be the same subject id in the genotype data for that subject
   to be included in tha analysis. 

  \bold{Missing data:} If any of the variables defined in \code{cc.var} or \code{covars.obj}
   contain missing values, then those subjects
  will be removed from the covariate and outcome data.  
}
\keyword{misc}


