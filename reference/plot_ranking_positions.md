# Plot how often items appear at each ranking position

This descriptive heatmap is useful before fitting a PL model. It shows
which items commonly appear near the top of observed rankings. It is not
adjusted for uncertainty or latent grouping.

## Usage

``` r
plot_ranking_positions(rankings)
```

## Arguments

- rankings:

  A btsbm_ranking_data object created by as_rankings().

## Value

A ggplot object.
