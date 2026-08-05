# Compute pointwise log likelihoods from a BTSBM fit

Compute pointwise log likelihoods from a BTSBM fit

## Usage

``` r
log_lik(object, ...)
```

## Arguments

- object:

  A `btsbm_fit` object.

- ...:

  Unused.

## Value

An `S × R` matrix with one column for each pair-count cell (BT) or whole
ranking row (PL). It has a `unit_index` attribute.
