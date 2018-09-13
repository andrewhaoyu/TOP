
\name{modelInfo}
\alias{modelInfo}
\title{Model Information}
\description{
 Returns objects based on the model information and options
}
\usage{
modelInfo(data, cc.var, subtype.vars, covars.obj=NULL, options=NULL)
}
\arguments{
  \item{data}{A data frame. No default.}
  \item{cc.var}{Name of the numeric case-control column.
    This column must be coded as 0-1, where 0 denotes a 
    control. Subjects with missing values in this column 
   will be removed from the analysis. No default.}
  \item{subtype.vars}{Character vector of column names for
   the subtype variables. This vector must be of length
  greater than one. See details. No default.}
  \item{covars.obj}{A formula or character vector of variable names 
     that define the covariate adjustments. This object should not 
     contain any of \code{test.vars}.
     An example formula is \code{~gender + age + ev1 + ev2}.
     No default.}
  \item{options}{List of options. See \code{\link{options}}}.
}
\value{
   A list with names \code{stage1to2.model}, \code{outcomes},
  \code{outcomes.removed}, \code{delta0}, and \code{options}.
  The matrix \code{stage1to2.model} is the matrix that was used to link the first stage
  parameters to the second stage parameters.
  The matrix \code{outcomes} gives the unique outcomes defined by \code{subtype.vars}
  and the integer value associated with each.
  The object \code{outcomes.removed} will be a matrix of the outcomes removed from
  the analysis or \code{NULL} if none of them were removed.
  The vector \code{delta0} contains the initial parameter estimates.
}
\details{
  This function can be called before an analysis is done to see if any outcomes
  will be removed due to small frequency counts, or to see the matrix that maps
  the first stage parameters ro the second stage parameters.
}

%\references{}

\seealso{ \code{\link{fixedEffect}}, \code{\link{randomEffect}}}

\examples{
  data(data, package="TOP")
  modelInfo(data, "Behavior", c("ER", "PR"), covars.obj=~PC1 + PC2)
 }
 

\keyword{ regression, two-stage } 
