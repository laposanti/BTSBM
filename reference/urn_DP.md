# Urn weight: Dirichlet Process

Urn weight: Dirichlet Process

## Usage

``` r
urn_DP(v_minus, alpha_PY)
```

## Arguments

- v_minus:

  integer sizes of occupied clusters (excluding focal item).

- alpha_PY:

  numeric concentration (\>0).

## Value

Numeric vector length \\K+1\\: existing weights and new-cluster mass.
\#' @references Legramanti, S., Rigon, T., Durante, D., Dunson, D.B.,
2022. Extended stochastic block models with application to criminal
networks. The Annals of Applied Statistics 16.
https://doi.org/10.1214/21-aoas1595
