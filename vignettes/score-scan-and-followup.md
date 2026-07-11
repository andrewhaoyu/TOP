<!-- This file is generated from `score-scan-and-followup.Rmd`. Edit the source and re-render it; do not edit this generated file directly. -->

Score scans and detailed follow-up
================

For a large variant scan, fit the nuisance model once, calculate score
tests in batches, and reserve full TOP models for selected variants.

``` r
library(TOP)
data("data", package = "TOP")

y <- data[, 1:5]
variants <- data[, 6, drop = FALSE]
adjustment <- data[, 7, drop = FALSE]
```

## Fit reusable score support

`top_score_support()` fits the null model containing adjustment
covariates. The returned object can be reused for every batch with the
same subjects, outcomes, and nuisance model.

``` r
support <- top_score_support(
  y = y,
  additive = adjustment,
  missing_value = 888
)

support
#> TOP score-scan support
#> Observations: 10000  Retained subtypes: 21 
#> Dropped observations: 11
```

## Scan variants

Each variant column is tested independently. The default additive
structure tests the baseline effect and all tumor-marker main effects
jointly.

``` r
scan <- top_score_scan(
  y = y,
  variants = variants,
  support = support,
  structure = "additive",
  missing_value = 888
)

as.data.frame(scan)
#>   variant statistic df    p.value status
#> 1     SNP  10.91802  5 0.05302981     ok
scan$score
#>     baseline effect ER main effect PR main effect HER2 main effect
#> SNP        57.47495       46.01162        37.4334         53.06665
#>     grade main effect
#> SNP          107.4106
```

For genome-wide data, pass manageable column batches while reusing
`support`.

``` r
all_results <- lapply(variant_batches, function(batch) {
  as.data.frame(top_score_scan(
    y = y,
    variants = batch,
    support = support,
    missing_value = 888
  ))
})
all_results <- do.call(rbind, all_results)
```

Singular information matrices are recorded in the `status` column with
an `NA` p-value, allowing a batch scan to continue.

## Fit selected variants

`top_followup()` fits one detailed Wald model per selected variant and
extracts variant-specific global association, global heterogeneity, and
second-stage effect results.

``` r
followup <- top_followup(
  y = y,
  variants = variants,
  adjustment = adjustment,
  variant_structure = "additive",
  missing_value = 888,
  keep_models = FALSE
)

followup$results
#>   variant global_association_p global_heterogeneity_p status
#> 1     SNP               0.0543                  0.126     ok
followup$effects
#>   Covariate SecondStageEffect OddsRatio OddsRatio(95%CI low)
#> 1       SNP   baseline effect      0.90                 0.76
#> 2       SNP    ER main effect      0.95                 0.85
#> 3       SNP    PR main effect      1.05                 0.93
#> 4       SNP  HER2 main effect      1.16                 1.01
#> 5       SNP grade main effect      1.08                 1.01
#>   OddsRatio(95%CI high) Pvalue
#> 1                  1.07 0.2390
#> 2                  1.08 0.4550
#> 3                  1.19 0.4490
#> 4                  1.35 0.0400
#> 5                  1.17 0.0368
```

Use `on_error = "record"` (the default) for several selected variants so
a single nonconvergent model is reported without discarding successful
fits. Use `keep_models = TRUE` when full fitted `top_model` objects are
needed.
