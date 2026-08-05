# Compute implied item abilities from a fit

`implied_ability()` is retained as a compatibility alias for
[`posterior_strength()`](https://laposanti.github.io/BTSBM/reference/posterior_strength.md).
New code should use
[`posterior_strength()`](https://laposanti.github.io/BTSBM/reference/posterior_strength.md).

## Usage

``` r
implied_ability(fit, summary = c("draws", "mean"))
```

## Arguments

- fit:

  A simple BT/PL or item-SBM fit.

- summary:

  Either `"draws"` or `"mean"`.

## Value

A draws matrix or named posterior mean vector.
