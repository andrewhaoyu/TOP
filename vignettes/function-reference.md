# TOP function reference

This is a practical map of the supported TOP interface. In R, open the full
argument-level help with `?function_name`, for example `?top_fit`.

## Recommended workflow

| Function | Use it for |
| --- | --- |
| `top_fit()` | Fit one validated two-stage polytomous logistic-regression model. |
| `top_fit_custom()` | Fit a model with a scientifically defined, custom second-stage design. |
| `top_score_support()` | Fit and retain the null-model support used for a multi-variant score scan. |
| `top_score_scan()` | Score-test a matrix of variants, reusing `top_score_support()` output. |
| `top_followup()` | Fit detailed models for selected variants after a scan. |
| `top_meta()` | Combine aligned score vectors and information matrices from independent studies. |

`top_fit()` returns a `top_model`. Its principal components are
`second.stage.test`, `first.stage.test`, and `global.test`; `coef()`,
`vcov()`, `summary()`, and `print()` are available for the fitted object.

## Input and design helpers

| Function | Use it for |
| --- | --- |
| `GenerateZstandard()` | Create TOP's standard second-stage design from a tumor-characteristic matrix. |
| `GenerateZDesignBaselineonly()` | Construct a baseline-only design. |
| `GenerateZDesignAdditive()` | Construct a design with marker main effects. |
| `GenerateZDesignPairwiseInteraction()` | Construct a design with marker main effects and pairwise interactions. |
| `GenerateZDesignSaturated()` | Construct a design with a separate effect for each retained subtype. |
| `GenerateZDesignCombination()` | Construct a specified combination of second-stage effects. |

For the standard `top_*()` functions, supply an outcome `y` whose first
column is case-control status (`0`/`1`) and whose remaining columns are tumor
characteristics. Controls must have `NA` tumor characteristics. Supply a
numeric `missing_value` when cases use a distinct missing-tumor code.

## Legacy API

Existing scripts can continue to use `TwoStageModel()`, `ScoreTest()`,
`ScoreTestMixedModel()`, and `ScoreMetaAnalysis()`. New work should prefer the
`top_*()` functions above: they validate inputs and return consistently named
objects. The historical low-level helpers remain available for reproducibility,
but they are not the recommended starting point for new analyses.

For the complete legacy reference, browse the package's
[Rd source files](../man/).
