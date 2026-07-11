# Combine TOP score statistics across studies

Combines independent study-level score vectors and information matrices.
Any number of studies may be supplied, and study order does not affect
the result.

## Usage

``` r
top_meta(studies)
```

## Arguments

- studies:

  A non-empty list. Each element must be a list containing a numeric
  \`score\` vector and a square numeric \`information\` matrix.

## Value

An object of class \`top_meta\` containing the combined score,
information matrix, chi-squared statistic, degrees of freedom, p-value,
and number of studies.

## Examples

``` r
studies <- list(
  list(score = c(1, 2), information = diag(c(2, 3))),
  list(score = c(0.5, -0.5), information = diag(c(1, 2)))
)
top_meta(studies)
#> TOP score meta-analysis
#> Studies: 2 
#> Chi-squared: 1.2 on 2 df; p = 0.5488 
```
