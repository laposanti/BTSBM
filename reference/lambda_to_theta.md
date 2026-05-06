# Map \\\lambda\\ to Bradley–Terry \\\theta = \lambda_i / (\lambda_i + \lambda_j)\\

Map \\\lambda\\ to Bradley–Terry \\\theta = \lambda_i / (\lambda_i +
\lambda_j)\\

## Usage

``` r
lambda_to_theta(lambda)
```

## Arguments

- lambda:

  Numeric vector of positive rates \\\lambda_i\\.

## Value

A matrix with entries \\\theta\_{ij} = \lambda_i / (\lambda_i +
\lambda_j)\\. Diagonal is 1/2 by the formula; you may overwrite if you
prefer NA on the diagonal.

## Examples

``` r
lambda_to_theta(c(1,2,4))
#>           [,1]      [,2]      [,3]
#> [1,] 0.5000000 0.3333333 0.2000000
#> [2,] 0.6666667 0.5000000 0.3333333
#> [3,] 0.8000000 0.6666667 0.5000000
```
