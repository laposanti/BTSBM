# Urn weight: Dirichlet–Multinomial with finite cap K_DM

Urn weight: Dirichlet–Multinomial with finite cap K_DM

## Usage

``` r
urn_DM(v_minus, beta_DM, K_DM)
```

## Arguments

- v_minus:

  integer sizes of occupied clusters (excluding focal item).

- beta_DM:

  numeric concentration (\>0).

- K_DM:

  integer maximum number of clusters (\>=1).

## Value

Numeric vector length \\K+1\\: existing weights and new-cluster mass.
\#' @references Legramanti, S., Rigon, T., Durante, D., Dunson, D.B.,
2022. Extended stochastic block models with application to criminal
networks. The Annals of Applied Statistics 16.
https://doi.org/10.1214/21-aoas1595
