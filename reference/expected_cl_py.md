# Expected number of clusters under finite/inf PY/DM-like priors (helper)

Expected number of clusters under finite/inf PY/DM-like priors (helper)

## Usage

``` r
expected_cl_py(n, sigma, theta, K_DM)
```

## Arguments

- n:

  integer(1) number of items.

- sigma:

  numeric(1) discount (0 for DP/DM).

- theta:

  numeric(1) concentration parameter.

- K_DM:

  integer(1) maximum number of clusters, or `Inf`.

## Value

Numeric(1) expected number of clusters. \#' @references Legramanti, S.,
Rigon, T., Durante, D., Dunson, D.B., 2022. Extended stochastic block
models with application to criminal networks. The Annals of Applied
Statistics 16. https://doi.org/10.1214/21-aoas1595
