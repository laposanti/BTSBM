# Fit a Bayesian Bradley–Terry or Plackett–Luce model

`fit_btsbm()` is the common fitting entry point. This first
implementation supports simple BT, BT–SBM, simple PL, PL–SBM, PL ranking
mixtures, and PL–LBM. The latter two start with readable R Gibbs kernels
derived from the PLuce exponential-race architecture; their
high-performance C++ kernels are a later parity-tested optimisation.

## Usage

``` r
fit_btsbm(data, model, control = mcmc_control())
```

## Arguments

- data:

  A `btsbm_pairwise_data` or `btsbm_ranking_data` object created by
  [`as_bt_data()`](https://laposanti.github.io/BTSBM/reference/as_bt_data.md)
  or
  [`as_rankings()`](https://laposanti.github.io/BTSBM/reference/as_rankings.md).

- model:

  A `btsbm_model` object.

- control:

  A `btsbm_mcmc_control` object.

## Value

An object of class `btsbm_fit`.

## References

Caron, F. and Doucet, A. (2012). Efficient Bayesian inference for
generalized Bradley–Terry models. *Journal of Computational and
Graphical Statistics*, 21(1), 174–196.
[doi:10.1080/10618600.2012.638220](https://doi.org/10.1080/10618600.2012.638220)
