# Scan variants with a TOP fixed-effect score test

Tests each variant independently using a null model that is fitted once.
For large scans, create \`support\` with \[top_score_support()\] and
reuse it across variant batches.

## Usage

``` r
top_score_scan(
  y,
  variants,
  structure = c("additive", "baseline", "pairwise", "saturated"),
  support = NULL,
  baseline = NULL,
  additive = NULL,
  pairwise = NULL,
  saturated = NULL,
  missing_value = NULL,
  subtype_cutoff = 10
)
```

## Arguments

- y:

  Numeric outcome matrix as described in \[top_fit()\].

- variants:

  Numeric vector, matrix, or data frame with one variant per column and
  one subject per row.

- structure:

  Second-stage structure used for the tested variants.

- support:

  Optional object returned by \[top_score_support()\].

- baseline, additive, pairwise, saturated:

  Nuisance covariates used only when \`support\` is not supplied.

- missing_value:

  Numeric missing-tumor code for cases, or \`NULL\`.

- subtype_cutoff:

  Subtypes with counts at or below this value are omitted when fitting a
  new support object.

## Value

A \`top_score_scan\` object containing a tidy result table, score
matrix, information array, and reusable support object.
