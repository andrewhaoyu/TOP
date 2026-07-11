<!-- This file is generated from `meta-analysis.Rmd`. Edit the source and re-render it; do not edit this generated file directly. -->

Score-statistic meta-analysis
================

`top_meta()` combines score vectors and information matrices from any
number of independent studies. This supports distributed analyses in
which individual participant data cannot be pooled.

## Required study-level objects

Each study supplies a list with:

- `score`: a numeric score vector;
- `information`: the corresponding square symmetric information matrix.

All studies must use the same parameter order and second-stage design.

``` r
library(TOP)

studies <- list(
  array_a = list(
    score = c(1.0, 2.0),
    information = diag(c(2, 3))
  ),
  array_b = list(
    score = c(0.5, -0.5),
    information = diag(c(1, 2))
  ),
  array_c = list(
    score = c(-0.25, 1.0),
    information = diag(c(4, 1))
  )
)

meta <- top_meta(studies)
meta
#> TOP score meta-analysis
#> Studies: 3 
#> Chi-squared: 1.2649 on 2 df; p = 0.5313
```

``` r
meta$score
#> [1] 1.25 2.50
meta$information
#>      [,1] [,2]
#> [1,]    7    0
#> [2,]    0    6
meta$statistic
#> [1] 1.264881
meta$p.value
#> [1] 0.5312936
```

The combined score and information are sums across studies. TOP then
evaluates the quadratic score statistic using the combined information
matrix.

## Defensive checks

`top_meta()` rejects non-finite scores, non-square or asymmetric
information matrices, inconsistent parameter dimensions, and a singular
combined information matrix. Align parameter names and allele coding
before constructing the list; study-specific column mapping belongs in
the analysis project rather than in the package.

`ScoreMetaAnalysis()` remains available for older two-study scripts, but
new code should use `top_meta()`.
