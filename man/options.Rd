% --- Source file: options.Rd ---
\name{options}
\alias{options}
\docType{data}
\title{
 List of options
}
\description{
 The list of options for two-stage ordinal logistic regression.
}
\format{
      The format is: List of 8
       \describe{
           \item{stage1to2.model}{A matrix or one of "additive", "baselineonly", "pairwise.interaction", 
             "saturated" to describe the link between the first stage and second stage parameters
             (see details).  
             The default is "additive".}
           \item{tol}{Stopping tolerance for the EM algorithm. The default is 1e-4.}
           \item{tolMaxstep}{Stopping tolerance for the maximization step of the EM algorithm. 
                 The default is 1e-4.}
           \item{maxiter}{Maximum number of iterations. The default is 100.}
           \item{min.ncase}{The minimum number of complete cases for each outcome. 
               The default is 10.}
           \item{min.ncontrol}{The minimum number of controls. The default is 10.}
           \item{delta0}{NULL or a vector of initial estimates for the second stage parameters. 
           See details for the order of this vector. The default is NULL.}
           \item{test.cols}{An integer valued vector giving the indices for the parameters
                 of interest to test.
             This option is only valid when \code{stage1to2.model} is a matrix, and must
                 be specified. See details.}
     }
 }
       
\details{ Let there be \code{M} different outcomes labeled 1, ..., M and 
    \code{N} covariates labeled \code{Cov1}, ..., \code{CovN}.
    The matrix that links the first stage parameters to the second stage
   parameters has the form:

  \tabular{rrrrr}{
                \tab delta1 \tab delta2 \tab ... \tab deltaK \cr
  (Intercept.1) \tab x      \tab x      \tab ... \tab x      \cr
  (Test.1)      \tab x      \tab x      \tab ... \tab x      \cr
  Cov1.1        \tab x      \tab x      \tab ... \tab x      \cr
  Cov2.1        \tab x      \tab x      \tab ... \tab x      \cr
  ...           \tab ...    \tab ...    \tab ... \tab ...    \cr
  CovN.1        \tab x      \tab x      \tab ... \tab x      \cr
  (Intercept.2) \tab x      \tab x      \tab ... \tab x      \cr
  (Test.2)      \tab x      \tab x      \tab ... \tab x      \cr
  Cov1.2        \tab x      \tab x      \tab ... \tab x      \cr
  Cov2.2        \tab x      \tab x      \tab ... \tab x      \cr
  ...           \tab ...    \tab ...    \tab ... \tab ...    \cr
  CovN.2        \tab x      \tab x      \tab ... \tab x      \cr
  ...           \tab ...    \tab ...    \tab ... \tab ...    \cr
  (Intercept.M) \tab x      \tab x      \tab ... \tab x      \cr
  (Test.M)      \tab x      \tab x      \tab ... \tab x      \cr
  Cov1.M        \tab x      \tab x      \tab ... \tab x      \cr
  Cov2.M        \tab x      \tab x      \tab ... \tab x      \cr
  ...           \tab ...    \tab ...    \tab ... \tab ...    \cr
  CovN.M        \tab x      \tab x      \tab ... \tab x      
  } % END: tabular
 
  The rownames of the above matrix are the first stage parameters, and the second
  stage parameters are the column names. The names \code{(Intercept.j)} are the
  intercepts for each outcome, and the names \code{(Test.j)} are the first
  stage parameters of interest among each outcome for each variable to test.
  Currently, if the user does not set \code{stage1to2.model} to a matrix, then
  there is no modeling performed on the intercept parameters. That is, there 
  will be \code{M} second stage intercept parameters. 

  If the user inputs a matrix for the option
  \code{stage1to2.model}, then it must have the same form as above, and the
  user must specify \code{test.cols}.
  The order of the vector \code{delta0} has the same order as the columns of
  the matrix above. The option \code{test.cols} specifies the column numbers  
  of the above matrix that are for the parameters of interest to test. 
 
}
%\examples{}
\keyword{misc}

