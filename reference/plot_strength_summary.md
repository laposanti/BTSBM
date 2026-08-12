# Plot posterior item strengths and their uncertainty

This forest plot answers a first analysis question: which items appear
stronger, and how uncertain is that ordering? Overlapping intervals mean
that the data do not cleanly distinguish the corresponding strengths.

## Usage

``` r
plot_strength_summary(fit, credible_mass = 0.9)
```

## Arguments

- fit:

  A simple BT/PL or item-SBM btsbm_fit object.

- credible_mass:

  Width of the central posterior interval, in (0, 1).

## Value

A ggplot object.
