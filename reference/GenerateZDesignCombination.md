# Generate the three Z design matrix given the tumor characteristics

Generate the three Z design matrix given the tumor characteristics

## Usage

``` r
GenerateZDesignCombination(y, missingTumorIndicator = 888, cutoff = 10)
```

## Arguments

- y:

  the phenotype file. The first column is the case control disease
  status. The other columns are the tumor characteristics status

- missingTumorIndicator:

  Numeric code for missing case tumor markers.

- cutoff:

  Subtypes with counts at or below this value are removed.

## Value

A result from this computation.
