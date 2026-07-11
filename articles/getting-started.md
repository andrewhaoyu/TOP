# Getting started with TOP

TOP fits two-stage polytomous logistic regression models for
case-control studies in which cases are characterized by multiple
disease or tumor markers. The recommended interface for new analyses is
[`top_fit()`](https://andrewhaoyu.github.io/TOP/reference/top_fit.md).

## Data layout

The outcome matrix has case-control status in its first column (`0` for
controls and `1` for cases) and tumor characteristics in the remaining
columns. Tumor characteristics must be `NA` for controls. If case tumor
characteristics can be missing, use a distinct numeric code and pass it
as `missing_value`.

``` r

library(TOP)
data("data", package = "TOP")

y <- data[, 1:5]
snp <- data[, 6, drop = FALSE]
pc1 <- data[, 7, drop = FALSE]

head(y)
#>      case.control ER PR HER2 grade
#> [1,]            0 NA NA   NA    NA
#> [2,]            0 NA NA   NA    NA
#> [3,]            1  1  0    1     1
#> [4,]            1  1  1  888   888
#> [5,]            1  1  1    0     0
#> [6,]            0 NA NA   NA    NA
table(y[, 1])
#> 
#>    0    1 
#> 4888 5112
```

The simulated package data use `888` for missing case tumor
characteristics.

## Fit an additive model

Here the SNP and principal component both use an additive second-stage
structure.
[`top_fit()`](https://andrewhaoyu.github.io/TOP/reference/top_fit.md)
validates dimensions, coding, missingness, and finite covariate values
before calling the numerical model.

``` r

fit <- top_fit(
  y = y,
  additive = cbind(snp, pc1),
  missing_value = 888
)

fit
#> TOP two-stage polytomous model
#> Observations: 10000  Tumor markers: 4 
#> Log-likelihood: -16726.5  AIC: 33514.9
```

Results have stable names. The second-stage table contains baseline
case-control effects and marker-specific case-case effects.

``` r

fit$second.stage.test
#>    Covariate SecondStageEffect OddsRatio OddsRatio(95%CI low)
#> 1        SNP   baseline effect      0.90                 0.76
#> 2        SNP    ER main effect      0.95                 0.85
#> 3        SNP    PR main effect      1.05                 0.93
#> 4        SNP  HER2 main effect      1.16                 1.01
#> 5        SNP grade main effect      1.08                 1.01
#> 6        PC1   baseline effect      1.10                 0.99
#> 7        PC1    ER main effect      0.99                 0.92
#> 8        PC1    PR main effect      0.95                 0.88
#> 9        PC1  HER2 main effect      0.95                 0.87
#> 10       PC1 grade main effect      1.03                 0.99
#>    OddsRatio(95%CI high) Pvalue
#> 1                   1.07 0.2390
#> 2                   1.08 0.4550
#> 3                   1.19 0.4490
#> 4                   1.35 0.0400
#> 5                   1.17 0.0368
#> 6                   1.22 0.0688
#> 7                   1.07 0.8300
#> 8                   1.02 0.1720
#> 9                   1.04 0.2340
#> 10                  1.08 0.1800
fit$global.test
#>   Covariate Global test for association Global test for heterogeneity
#> 1       SNP                     0.05430                         0.126
#> 2       PC1                     0.00563                         0.251
head(fit$first.stage.test)
#>   Covariate          Subtypes OddsRatio OddsRatio(95%CI low)
#> 1       SNP ER1PR0HER20grade0      0.86                 0.71
#> 2       PC1 ER1PR0HER20grade0      1.09                 0.97
#> 3       SNP ER0PR1HER20grade0      0.95                 0.84
#> 4       PC1 ER0PR1HER20grade0      1.04                 0.97
#> 5       SNP ER1PR1HER20grade0      0.90                 0.77
#> 6       PC1 ER1PR1HER20grade0      1.04                 0.94
#>   OddsRatio(95%CI high) Pvalue
#> 1                  1.05  0.137
#> 2                  1.23  0.143
#> 3                  1.07  0.377
#> 4                  1.12  0.252
#> 5                  1.06  0.223
#> 6                  1.14  0.487
```

Standard R methods are also available.

``` r

head(coef(fit, stage = "second"))
#> [1] -3.581787 -3.799975 -4.021544 -4.638728 -3.177350 -4.688342
dim(vcov(fit, stage = "second"))
#> [1] 31 31
summary(fit)
#> TOP model summary
#> Call: top_fit(y = y, additive = cbind(snp, pc1), missing_value = 888)
#> Log-likelihood: -16726.5  AIC: 33514.9 
#> 
#> Global tests:
#>   Covariate Global test for association Global test for heterogeneity
#> 1       SNP                     0.05430                         0.126
#> 2       PC1                     0.00563                         0.251
```

## Choose a second-stage structure

Covariates can be assigned to different structures:

- `baseline`: a common association across subtypes;
- `additive`: baseline plus tumor-marker main effects;
- `pairwise`: additive effects plus pairwise marker interactions;
- `saturated`: a separate effect for every retained subtype.

For example, a common SNP effect can be adjusted for an additive PC
effect.

``` r

fit_mixed <- top_fit(
  y = y,
  baseline = snp,
  additive = pc1,
  missing_value = 888
)
```

## Subtype counts and convergence

TOP historically omits subtypes with counts at or below the cutoff. The
default `subtype_cutoff = 10` preserves that behavior. Lower cutoffs can
make information matrices unstable, particularly in small datasets.
Numerical failures raise R errors instead of terminating the R session.

Set `options(TOP.verbose = TRUE)` temporarily when EM progress messages
are useful.

## Existing scripts

[`TwoStageModel()`](https://andrewhaoyu.github.io/TOP/reference/TwoStageModel.md)
remains available for backward compatibility. New code should prefer
[`top_fit()`](https://andrewhaoyu.github.io/TOP/reference/top_fit.md)
and named components such as `fit$second.stage.test` rather than
positional indexing.
