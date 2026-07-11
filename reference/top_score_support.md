# Prepare reusable null-model support for TOP score scans

Fits the nuisance model once and stores the matrices needed to scan many
variants. Reuse the returned object in \[top_score_scan()\] for
successive variant batches.

## Usage

``` r
top_score_support(
  y,
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

- baseline, additive, pairwise, saturated:

  Numeric nuisance covariates for each second-stage structure.

- missing_value:

  Numeric missing-tumor code for cases, or \`NULL\`.

- subtype_cutoff:

  Subtypes with counts at or below this value are omitted.

## Value

A reusable object of class \`top_score_support\`.
