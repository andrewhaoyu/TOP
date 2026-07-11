# Fit a TOP model with validated inputs and named results

This is the recommended interface for new analyses. The original
\[TwoStageModel()\] function remains available for backward
compatibility.

## Usage

``` r
top_fit(
  y,
  baseline = NULL,
  additive = NULL,
  pairwise = NULL,
  saturated = NULL,
  missing_value = NULL,
  start = NULL,
  subtype_cutoff = 10
)
```

## Arguments

- y:

  Numeric matrix or data frame whose first column is case-control status
  (0/1) and remaining columns are tumor characteristics. Tumor
  characteristics must be \`NA\` for controls.

- baseline:

  Covariates with a common baseline effect across subtypes.

- additive:

  Covariates using the additive second-stage model.

- pairwise:

  Covariates using the pairwise-interaction model.

- saturated:

  Covariates using the saturated model.

- missing_value:

  Numeric code used for missing tumor characteristics in cases, or
  \`NULL\` when case tumor data are complete.

- start:

  Optional starting vector for second-stage parameters.

- subtype_cutoff:

  Subtypes with counts at or below this value are omitted, matching the
  historical TOP behavior.

## Value

A named list of class \`top_model\`.
