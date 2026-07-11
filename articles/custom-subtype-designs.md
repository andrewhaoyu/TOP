# Custom subtype designs

[`top_fit_custom()`](https://andrewhaoyu.github.io/TOP/reference/top_fit_custom.md)
supports scientifically defined subtype groupings that are not
represented by TOP’s predefined hierarchical structures.

``` r

library(TOP)
data("data", package = "TOP")

y <- data[, 1:5]
snp <- data[, 6, drop = FALSE]
pc1 <- data[, 7, drop = FALSE]
z_standard <- GenerateZstandard(y, missingTumorIndicator = 888)
```

Rows of the custom design correspond to the retained first-stage
subtypes. Columns define custom second-stage effects. This example maps
the cross-classified subtypes into five groups based on hormone-receptor
status, HER2, and grade.

``` r

design <- matrix(0, nrow(z_standard), 5)
colnames(design) <- c(
  "HR+_HER2-_lowgrade",
  "HR+_HER2+",
  "HR+_HER2-_highgrade",
  "HR-_HER2+",
  "HR-_HER2-"
)

hr_positive <- z_standard[, 1] == 1 | z_standard[, 2] == 1
her2_positive <- z_standard[, 3] == 1
high_grade <- z_standard[, 4] == 2

design[hr_positive & !her2_positive & !high_grade, 1] <- 1
design[hr_positive & her2_positive, 2] <- 1
design[hr_positive & !her2_positive & high_grade, 3] <- 1
design[!hr_positive & her2_positive, 4] <- 1
design[!hr_positive & !her2_positive, 5] <- 1

table(rowSums(design))
#> 
#>  1 
#> 21
```

Every retained subtype must map to at least one custom effect.
Overlapping columns are allowed when scientifically intended, provided
the resulting model is identifiable.

``` r

custom_fit <- top_fit_custom(
  y = y,
  covariate = snp,
  design = design,
  adjustment = pc1,
  missing_value = 888
)

custom_fit$second.stage.test
#>    Covariate   SecondStageEffect OddsRatio OddsRatio(95%CI low)
#> 1        SNP  HR+_HER2-_lowgrade      1.02                 0.93
#> 2        SNP           HR+_HER2+      1.12                 1.02
#> 3        SNP HR+_HER2-_highgrade      1.06                 0.93
#> 4        SNP           HR-_HER2+      1.03                 0.73
#> 5        SNP           HR-_HER2-      0.94                 0.66
#> 6        PC1     baseline effect      1.10                 0.99
#> 7        PC1      ER main effect      0.99                 0.92
#> 8        PC1      PR main effect      0.95                 0.88
#> 9        PC1    HER2 main effect      0.95                 0.87
#> 10       PC1   grade main effect      1.03                 0.99
#>    OddsRatio(95%CI high) Pvalue
#> 1                   1.12 0.7050
#> 2                   1.22 0.0158
#> 3                   1.22 0.3860
#> 4                   1.44 0.8740
#> 5                   1.36 0.7580
#> 6                   1.22 0.0710
#> 7                   1.07 0.8050
#> 8                   1.02 0.1720
#> 9                   1.04 0.2450
#> 10                  1.08 0.1690
custom_fit$global.test
#>   Covariate Global test for association Global test for heterogeneity
#> 1       SNP                     0.24800                         0.542
#> 2       PC1                     0.00548                         0.246
```

The fitted object uses the same `top_model` methods as
[`top_fit()`](https://andrewhaoyu.github.io/TOP/reference/top_fit.md).

``` r

custom_fit
#> TOP two-stage polytomous model
#> Observations: 10000  Tumor markers: 4 
#> Log-likelihood: -16728.5  AIC: 33519
head(coef(custom_fit))
#> [1] -3.659895 -3.833142 -4.080667 -4.625928 -3.235701 -4.700390
```
