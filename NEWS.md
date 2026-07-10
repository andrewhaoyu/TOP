# TOP 1.1.0

## Build and native-code reliability

- Modernized 145 legacy C function definitions for C23-compatible compilers.
- Registered all native routines and disabled dynamic symbol lookup.
- Replaced native `exit(1)` calls with R errors and routed native output through
  the R console.

## User interface

- Added the validated `top_fit()` interface and `top_model` print, summary,
  coefficient, and covariance methods.
- Added `top_fit_custom()` for scientifically defined second-stage designs.
- Added reusable `top_score_support()`, batch-friendly `top_score_scan()`, and
  detailed `top_followup()` workflows.
- Added `top_meta()` for score-statistic meta-analysis across any number of
  studies; `ScoreMetaAnalysis()` remains backward compatible.
- Added stable names to the legacy model result components.
- Invalid control tumor coding now raises an error instead of returning text.

## Testing and documentation

- Replaced environment-specific analysis scripts with reproducible `testthat`
  tests while preserving the original scripts under `legacy/tests`.
- Added R CMD check CI for Linux, macOS, and Windows on devel, release, and
  old-release R.
- Added a pkgdown site workflow, four focused vignettes, GitHub issue forms,
  contribution guidance, citation metadata, a runnable README quick start,
  dataset documentation, and synchronized function documentation.
