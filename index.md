# TOP

[![R-CMD-check](https://github.com/andrewhaoyu/TOP/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/andrewhaoyu/TOP/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/andrewhaoyu/TOP/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/andrewhaoyu/TOP/actions/workflows/pkgdown.yaml)
[![License:
GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

TOP fits two-stage polytomous logistic regression models for
case-control studies with multiple disease or tumor characteristics. It
supports missing tumor-marker data, reusable fixed-effect score scans,
detailed variant follow-up, custom subtype designs, and score-statistic
meta-analysis.

## Installation

Install the development release from GitHub:

``` r

install.packages("remotes")
remotes::install_github("andrewhaoyu/TOP")
```

Then load TOP:

``` r

library(TOP)
```

## Quick start

The outcome matrix contains case-control status in column 1 and tumor
characteristics in the remaining columns. Controls must have `NA` tumor
characteristics. In the package example, `888` marks missing tumor data
among cases.

``` r

data("data", package = "TOP")

y <- data[, 1:5]
snp <- data[, 6, drop = FALSE]
pc1 <- data[, 7, drop = FALSE]

fit <- top_fit(
  y = y,
  additive = cbind(snp, pc1),
  missing_value = 888
)

fit
fit$second.stage.test
fit$global.test
summary(fit)
```

TOP model results have stable names and standard methods:

``` r

coef(fit, stage = "second")
vcov(fit, stage = "second")
```

## Genome-wide score scans

Fit the nuisance model once and reuse it across variant batches:

``` r

support <- top_score_support(
  y = y,
  additive = pc1,
  missing_value = 888
)

scan <- top_score_scan(
  y = y,
  variants = snp,
  support = support,
  missing_value = 888
)

as.data.frame(scan)
```

Fit detailed Wald models only for selected variants:

``` r

followup <- top_followup(
  y = y,
  variants = snp,
  adjustment = pc1,
  missing_value = 888,
  keep_models = FALSE
)

followup$results
followup$effects
```

## Meta-analysis

Combine aligned score vectors and information matrices from any number
of independent studies:

``` r

studies <- list(
  list(score = c(1, 2), information = diag(c(2, 3))),
  list(score = c(0.5, -0.5), information = diag(c(1, 2))),
  list(score = c(-0.25, 1), information = diag(c(4, 1)))
)

top_meta(studies)
```

## Documentation

- [Getting
  started](https://andrewhaoyu.github.io/TOP/vignettes/getting-started.md)
  ([source](https://andrewhaoyu.github.io/TOP/vignettes/getting-started.Rmd))
- [Score scans and
  follow-up](https://andrewhaoyu.github.io/TOP/vignettes/score-scan-and-followup.md)
  ([source](https://andrewhaoyu.github.io/TOP/vignettes/score-scan-and-followup.Rmd))
- [Score-statistic
  meta-analysis](https://andrewhaoyu.github.io/TOP/vignettes/meta-analysis.md)
  ([source](https://andrewhaoyu.github.io/TOP/vignettes/meta-analysis.Rmd))
- [Custom subtype
  designs](https://andrewhaoyu.github.io/TOP/vignettes/custom-subtype-designs.md)
  ([source](https://andrewhaoyu.github.io/TOP/vignettes/custom-subtype-designs.Rmd))
- [Function
  reference](https://andrewhaoyu.github.io/TOP/vignettes/function-reference.md)

The guides are rendered Markdown files that open readably on GitHub. The
linked `.Rmd` files are their reproducible source. A richer package
website is also ready to publish from the `gh-pages` branch once GitHub
Pages is enabled for this repository.

[`TwoStageModel()`](https://andrewhaoyu.github.io/TOP/reference/TwoStageModel.md)
and the original score-test functions remain available for existing
scripts. New analyses should use the validated `top_*()` interfaces. The
original long methodological tutorial is archived under
[`legacy/tutorial`](https://andrewhaoyu.github.io/TOP/legacy/tutorial).

## Subtype filtering and convergence

TOP preserves its historical rule that a cutoff of 10 omits subtypes
with 10 or fewer cases. Sparse subtypes can produce unstable or singular
information matrices. Numerical failures now raise informative R errors
instead of terminating the R session. Set `options(TOP.verbose = TRUE)`
to show EM progress.

## Citation

If TOP contributes to an analysis, cite:

> Zhang H, Zhao N, Ahearn TU, Wheeler W, García-Closas M, Chatterjee N.
> A mixed-model approach for powerful testing of genetic associations
> with cancer risk incorporating tumor characteristics. *Biostatistics*.
> 2021;22(4):772–788.
> [doi:10.1093/biostatistics/kxz065](https://doi.org/10.1093/biostatistics/kxz065)

From R, run `citation("TOP")` for a formatted entry.

## Contributing

Bug reports and feature requests are welcome through [GitHub
Issues](https://github.com/andrewhaoyu/TOP/issues). See
[`CONTRIBUTING.md`](https://andrewhaoyu.github.io/TOP/CONTRIBUTING.md)
before proposing statistical or API changes. Never attach
participant-level or sensitive study data to a public issue.
