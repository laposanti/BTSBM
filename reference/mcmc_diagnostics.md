# Summarise MCMC mixing diagnostics

Computes a lightweight initial-positive-sequence effective sample size,
Monte Carlo standard error, and lag-one autocorrelation for saved scalar
traces. PL fits include the total ranking log likelihood; block models
also include their occupied-cluster traces.

## Usage

``` r
mcmc_diagnostics(fit)
```

## Arguments

- fit:

  A `btsbm_fit` object.

## Value

A data frame with one row per available scalar trace.
