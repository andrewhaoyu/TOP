# Generate the combinations of all the tumor characteristics.

Generate the combinations of all the tumor characteristics.

## Usage

``` r
GenerateZstandard(y, missingTumorIndicator = 888, cutoff = 10)
```

## Arguments

- y:

  the phenotype file. The first column is the case control disease
  status. The other columns are the tumor characteristics status

- missingTumorIndicator:

  The indicators to show the tumor characteristics are missing. In the
  example, we put missing tumor characteristics as 888. Note, for all
  the controls subjects, they don't have tumor characteristics. So their
  tumor characteristics are put as NA instead of 888 to differentiate
  with cases missing tumor characteristics.

- cutoff:

  subtypes with counts at or below this value are removed

## Value

a matrix with all the combinations of tumor characteristics. All the
subtypes with less than 10 are removed by default.

## Examples

``` r
data(data, package="TOP") #load in the breast cancer example
#this is a simulated breast cancer example
#there are around 5000 breast cancer cases and 5000 controls, i.e. people without disease
data[1:5,]
#>      case.control ER PR HER2 grade SNP        PC1
#> [1,]            0 NA NA   NA    NA   0 -0.4382773
#> [2,]            0 NA NA   NA    NA   0  0.4723048
#> [3,]            1  1  0    1     1   1  0.3012679
#> [4,]            1  1  1  888   888   1 -0.2916645
#> [5,]            1  1  1    0     0   1  0.7650199
#four different tumor characteristics were included, ER (positive vs negative), PR (positive vs negative), HER2 (positive vs negative), grade (ordinal 1, 2, 3)
#the phenotype file
y <- data[,1:5]
#generate the combinations of all the subtypes
#by default, we remove all the subtypes with less than 10 cases
z.standard <- GenerateZstandard(y)
```
