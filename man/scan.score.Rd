
\name{scan.score}
\alias{scan.score}
\title{Scan function }
\description{
  Function to apply the two-stage model on an entire file of genotype data. 
}
\usage{ scan.score(geno.list, pheno.list, scan.options=NULL)
}
\arguments{
  \item{geno.list}{ List describing the genotype data. 
                    See \code{\link{geno.list}}. No default }
  \item{pheno.list}{ List describing the phenotype data.
                   See \code{\link{pheno.list}}. No default. }
  \item{scan.options}{ List of options (see details). The default is NULL.}
}
\details{
 This function can be called to run single SNP analyses on a subset
  of SNPs or an entire file of genotype data. 
  To use this function, the data must be stored in files as defined in
\code{\link{geno.list}} and \code{\link{pheno.list}}. 
Currently, the format of the genotype data must be the format from the IMPUTE2 software
 containing genotype probabilities.
 The subjects to be included in the analysis are the subjects with matching ids in the phenotype and
 genotype data. \cr

  \bold{scan.options list:} 
  \itemize{
    \item \code{options} See \code{\link{options}}.
    \item \code{out.file} File name to save the results.
       The last row of this file will have a row of "*" to signify
        normal completion.
       If NULL, then a default output file will be created in the working 
       directory as: \cr
      \code{paste(getwd(), "/", basename(geno.list$file), "_out_", a, "_", b, sep="")}
     where \code{a = geno.list$start.row} and \code{b = geno.list$stop.row}. 
    \item \code{method} "fixed" or "random". 
           The default is "fixed".
    \item \code{min.MAF} A value betwen 0 and 1 to only process SNPs with
         MAF greater than or equal to \code{min.MAF}.
                 The default is 0.
    \item \code{delta0.list} A list with names \code{save.N} and \code{max.r2} to
               store SNPs and estimated parameters in order obtain better initial
               estimates for a future SNP (see initial estimates below). 
               This list is only used with \code{method} = "random".    
                 The default is \code{list(save.N=50, max.r2=0.5)}.
  } % END: itemize 

  \bold{initial estimates:} \cr
    The random effects model takes significantly more time than the fixed 
    effect method since all second stage
    parameters need to be estimated with each SNP.  
    In order to speed up the random effects method, good initial estimates
    are needed. One way to accomplish this is with \code{delta0.list} and take
    advantage of the LD structure of SNPs in a genotype file. The idea to save 
   a set of independent snps and the corresponding parameter estimates in matrices
   and use a previously estimated set of parameters as the intial estimates for a future
  SNP in the genotype file that is in high LD with a saved SNP. A SNP is saved
  if its maximum R-squared value with previously saved SNPs is less than 
  \code{max.r2}. A maximum number of \code{save.N} SNPs and their corresponding
  estimated parameters will be saved. If a new SNP needs to be saved when there
  are already \code{save.N} SNPs saved, then the matrix of saved SNPs will be 
  overwritten starting with the first saved SNP. If the user does not want to use
  initial estimates created this way, then the user can set \code{save.N} to 0 or
  specify the initial estimates using \code{options$delta0}.
    

  
}
\value{
  This function returns NULL. The results are stored in \code{out.file}. The output file
  will contain columns specifying the SNP, position, reference allele, effect allele,
  effect allele frequency, p-values, and any error message.
  If \code{method} = "random", then \code{out.file} will contain a p-value for
  the fixed effects global test of association (GTA) and a p-value for
  random effects GTA.
 The last row of the file should be a row of "*" to signify normal completion.
}
%\references{ }
%\author{ }
\examples{
  gfile <- system.file("sampleData", "chr1_1.impute.txt", package="TOP")
  sfile <- system.file("sampleData", "subjects.impute.txt", package="TOP")
  pfile <- system.file("sampleData", "pheno.txt", package="TOP")

  # Define the list for the genotype data. 
 geno.list <- list(file=gfile, format="impute", subject.list=sfile)
 
 # Only process the first 5 SNPs in the file
 geno.list$start.vec <- 1
 geno.list$stop.vec <- 6

 # Define pheno.list
 pheno.list <- list(file=pfile, delimiter="\t", id.var="ID")
 
 # Define the variables in the model: Status ~ Age
 pheno.list$response.var <- "Status"
 pheno.list$main.vars <- "Age"

 # scan.score(snp.list, pheno.list)

}
\keyword{ GWAS } 
\seealso{ \code{\link{geno.list}}, \code{\link{pheno.list}} }
