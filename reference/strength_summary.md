# Summarise posterior item strengths with credible intervals

This is a compact numerical companion to plot_strength_summary(). The
reported strengths are relative: for Plackett–Luce fits, the default
reporting scale has geometric mean one; for Bradley–Terry fits, a larger
value means a higher chance of beating an item with a smaller value.

## Usage

``` r
strength_summary(fit, credible_mass = 0.9)
```

## Arguments

- fit:

  A simple BT/PL or item-SBM btsbm_fit object.

- credible_mass:

  Width of the central posterior interval, in (0, 1).

## Value

A data frame with one row per item, ordered from larger to smaller
posterior mean strength.
