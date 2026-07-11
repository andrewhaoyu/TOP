# Fit detailed TOP models for selected follow-up variants

Fits one TOP model per selected variant and returns variant-specific
global tests and second-stage effect estimates. This is intended for the
limited set of variants selected after a score scan, not for genome-wide
fitting.

## Usage

``` r
top_followup(
  y,
  variants,
  adjustment = NULL,
  variant_structure = c("additive", "baseline", "pairwise", "saturated"),
  missing_value = NULL,
  start = NULL,
  subtype_cutoff = 10,
  on_error = c("record", "stop"),
  keep_models = TRUE
)
```

## Arguments

- y:

  Numeric outcome matrix as described in \[top_fit()\].

- variants:

  Numeric vector, matrix, or data frame with selected variants in
  columns.

- adjustment:

  Optional numeric adjustment covariates. These use the additive
  second-stage structure.

- variant_structure:

  Second-stage structure for each tested variant.

- missing_value:

  Numeric missing-tumor code for cases, or \`NULL\`.

- start:

  Optional starting vector passed to \[top_fit()\].

- subtype_cutoff:

  Subtypes with counts at or below this value are omitted.

- on_error:

  Whether a failed fit should be recorded or stop the run.

- keep_models:

  Keep the fitted \`top_model\` objects in the result.

## Value

A \`top_followup\` object with global tests, effect estimates, optional
fitted models, and recorded errors.
