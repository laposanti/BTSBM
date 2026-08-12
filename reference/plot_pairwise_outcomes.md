# Plot observed pairwise win counts

This descriptive heatmap is useful before fitting a Bradley–Terry model.
Rows are winners and columns are losers; it displays the number of
observed wins, not a fitted probability.

## Usage

``` r
plot_pairwise_outcomes(pairwise)
```

## Arguments

- pairwise:

  A btsbm_pairwise_data object created by as_bt_data().

## Value

A ggplot object.
