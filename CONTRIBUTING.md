# Contributing to TOP

Contributions that improve statistical reliability, portability,
performance, documentation, or analyst usability are welcome.

## Before opening a pull request

1.  Open an issue for changes that alter statistical behavior or the
    public API.
2.  Work from a focused branch based on `master`.
3.  Add a reproducible test using simulated or package data.
4.  Update the relevant help page or vignette.
5.  Run `R CMD build .` followed by `R CMD check TOP_*.tar.gz`.

Please do not include participant-level or otherwise sensitive study
data.

## Compatibility

The `top_*()` functions are the recommended API. Existing exported
functions are retained for backward compatibility. A change to
estimates, tests, convergence, missing-data handling, or subtype
filtering must be called out explicitly and supported by regression
tests.

## Native code

Native routines must be registered, callable through the R API, and must
never call `exit()`, write directly to stdout/stderr, or use
unregistered symbols.

## Reporting results

Include the operating system, R version, TOP version, a minimal example,
and [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html) when
reporting a failure.
