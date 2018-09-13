
\name{TOP}
\alias{TOP}
\docType{package}
\title{
 TwO-stage Polytmous logistic regression
}
\description{
 A two-stage regression model for epidemiological studies with multivariate disease classification data.
}
\details{
Polytomous logistic regression is commonly used to analyze epidemiological data with disease subtype information. 
In this approach, effects of exposures on different disease subtypes are studied through separate exposure odds 
ratios comparing different case groups to the common control group. 
This package considers the situation where disease subtypes can be defined using multiple characteristics of a disease. 
For efficient analysis of such data, a two-stage modeling approach is used. 
At the first stage, a standard polytomous logistic regression model is considered for all possible distinct disease subtypes 
that can be defined by the cross-classification of the different disease characteristics. 
At the second stage, the exposure odds ratio parameters for the first-stage disease subtypes are further modeled in terms
 of the defining characteristics of the subtypes. When the total number of first-stage disease subtypes is small, 
standard maximum likelihood methods can be used for inference in the proposed model. 
For dealing with a large number of disease subtypes, a novel semiparametric pseudo-conditional-likelihood approach 
is applied that does not require any model assumption about the baseline probabilities for the different disease subtypes. 
}
\references{ 

Chatterjee, N. A Two-Stage Regression Model for Epidemiological Studies with 
Multivariate Disease Classification Data.
  Journal of the American Statistical Association, Vol. 99, No. 465 (Mar., 2004), pp.127-138.
}
\author{Nilanjan Chatterjee and Haoyu Zhang}
\keyword{package}
