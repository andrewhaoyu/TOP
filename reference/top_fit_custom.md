# Fit a TOP model with a custom second-stage design

Use this interface when predefined baseline, additive, pairwise, or
saturated structures do not represent the desired subtype grouping.

## Usage

``` r
top_fit_custom(
  y,
  covariate,
  design,
  adjustment = NULL,
  missing_value = NULL,
  start = NULL,
  subtype_cutoff = 10
)
```

## Arguments

- y:

  Numeric outcome matrix as described in \[top_fit()\].

- covariate:

  Numeric vector or one-column matrix whose effect is mapped by
  \`design\`.

- design:

  Numeric matrix with one row per retained subtype and one column per
  custom second-stage effect.

- adjustment:

  Optional additive adjustment covariates.

- missing_value:

  Numeric missing-tumor code for cases, or \`NULL\`.

- start:

  Optional starting vector.

- subtype_cutoff:

  Subtypes with counts at or below this value are omitted.

## Value

A named object of class \`top_model\`. Summary tables include both the
custom covariate and any jointly fitted adjustment covariates.
