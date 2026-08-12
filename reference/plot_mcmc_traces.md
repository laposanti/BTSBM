# Plot saved MCMC traces

Trace plots make it possible to see whether saved MCMC draws move around
a stable region. They are a first diagnostic, not a proof of
convergence.

## Usage

``` r
plot_mcmc_traces(fit, quantities = NULL)
```

## Arguments

- fit:

  A btsbm_fit object.

- quantities:

  Optional character vector of trace names to show.

## Value

A ggplot object.
