TOP
=======
Two-stage polytomous logistic regression (TOP): in this package, we developed a two-stage polytomous regression framework to handle cancer data with multivariate tumor characteristics. In the first stage, a standard polytomous model is used to specify for all subtypes defined by the cross-classification of different markers. In the second stage, the subtype-specific case-control odds ratios are specified using a more parsimonious model based on the case-control odds ratio for a baseline subtype, and the case-case parameters associated with tumor markers. Further, to reduce the degrees-of-freedom, we allow to specify case-case parameters for additional markers using a random-effect model.

Usage
=======
[The 'TOP' vignette](https://github.com/andrewhaoyu/TOP/blob/master/inst/TOP.pdf) will provide a good start point for using TOP package.


Development 
=======
This R package is developed by Haoyu Zhang and William Wheeler, and maintained by Haoyu Zhang <andrew.haoyu@gmail.com>.

Installation
=======
To install the development version of TOP, it's easiest to use the 'devtools' package.

install.packages("devtools")  
library(devtools)  
install_github("andrewhaoyu/TOP")

References
=======
Zhang, H., Zhao, N., Ahearn, T.U, Wheeler W., García-Closas, M., Chatterjee, N., A mixed-model approach for powerful testing of genetic associations with cancer risk incorporating tumor characteristics (Submitted)


