# Plot the posterior similarity of items or rankings

A dark cell means that two entities were often assigned to the same
latent group across posterior draws. The plot is label-invariant, which
makes it safer to interpret than the raw MCMC cluster labels.

## Usage

``` r
plot_posterior_similarity(
  fit,
  target = c("item", "ranking"),
  order_rows = TRUE
)
```

## Arguments

- fit:

  A btsbm_fit object containing an item or ranking partition.

- target:

  Whether to plot latent groups of "item" or "ranking".

- order_rows:

  Whether to cluster the rows and columns for display.

## Value

A ggplot object.
