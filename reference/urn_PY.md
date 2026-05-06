# Urn weight: Pitman–Yor process

Urn weight: Pitman–Yor process

## Usage

``` r
urn_PY(v_minus, alpha_PY, sigma_PY)
```

## Arguments

- v_minus:

  integer sizes of occupied clusters (excluding focal item).

- alpha_PY:

  numeric concentration (\> -sigma_PY).

- sigma_PY:

  numeric discount in \[0,1).

## Value

Numeric vector length \\K+1\\. \#' @references Legramanti, S., Rigon,
T., Durante, D., Dunson, D.B., 2022. Extended stochastic block models
with application to criminal networks. The Annals of Applied Statistics
16. https://doi.org/10.1214/21-aoas1595
