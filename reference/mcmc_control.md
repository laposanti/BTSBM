# Set MCMC controls for a BTSBM fit

Set MCMC controls for a BTSBM fit

## Usage

``` r
mcmc_control(
  iter = 2000,
  warmup = floor(iter/2),
  seed = NULL,
  store = c("minimal", "augmented"),
  init_item_cluster = NULL,
  init_ranking_cluster = NULL,
  verbose = FALSE,
  progress_every = 1000L,
  partition_moves = NULL
)
```

## Arguments

- iter:

  Total MCMC iterations.

- warmup:

  Number of warmup iterations.

- seed:

  Optional integer seed.

- store:

  Whether to store the minimal state or latent augmentations.

- init_item_cluster:

  Optional initial item partition.

- init_ranking_cluster:

  Optional initial ranking partition for PL mixture and PL–LBM models.

- verbose:

  Whether to report progress.

- progress_every:

  Progress frequency when `verbose` is `TRUE`.

- partition_moves:

  Optional result of
  [`partition_moves()`](https://laposanti.github.io/BTSBM/reference/partition_moves.md).

## Value

An object of class `btsbm_mcmc_control`.
