
\name{fixedEffect}
\alias{fixedEffect}
\title{Two stage ordinal logistic regression with fixed effects}
\description{
Two stage ordinal logistic regression with fixed effects
}
\usage{
fixedEffect(data, cc.var, subtype.vars, test.vars, covars.obj=NULL, options=NULL)
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
  \item{test.vars}{Character vector of variables to test individually. No default.}
  \item{covars.obj}{A formula or character vector of variable names 
     that define the covariate adjustments. This object should not 
     contain any of \code{test.vars}.
     An example formula is \code{~gender + age + ev1 + ev2}.
     No default.}
  \item{options}{List of options. See \code{\link{options}}}.
}
\value{
   A list with names \code{results}, \code{stage1to2.model}, \code{outcomes},
  and \code{outcomes.removed}.
  The matrix \code{stage1to2.model} is the matrix that was used to link the first stage
  parameters to the second stage parameters.
  The matrix \code{outcomes} gives the unique outcomes defined by \code{subtype.vars}
  and the integer value associated with each.
  The object \code{outcomes.removed} will be a matrix of the outcomes removed from
  the analysis or \code{NULL} if none of them were removed.

}
\details{
Subjects that contain a missing value for \code{cc.var}, for a covariate,
 or for \code{test.vars} will be removed from the analysis.
The columns defined by \code{subtype.vars} are allowed to have
missing values and the missing values will be imputed.
However there must be at least \code{options$min.ncase} complete cases for each outcome,
where complete cases are defined to be cases with all non-missing values 
for \code{subtype.vars}.
The set of all possible case outcomes is defined by the \code{subtype.vars} with 
all non-missing values. 

It is much more efficient to call this function with \code{length(test.vars) > 1} 
  than to make several calls with \code{length(test.vars}) = 1, 
   and is more efficient if none of the \code{test.vars} contain missing values.

}

%\references{}

\seealso{ \code{\link{randomEffect}}, \code{\link{options}} }

\examples{
  data(data, package="TOP")
  fixedEffect(data, "Behavior", c("ER", "PR"), "SNP", covars.obj=~PC1 + PC2)
 }
 

\keyword{ regression, two-stage } 
