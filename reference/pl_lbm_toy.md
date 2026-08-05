# Structured Plackett–Luce latent-block-model toy data

Generates rankings from a two-way latent block model with known
ranking-row clusters, item blocks, and a `C × K` latent-strength matrix.
It is a compact test and vignette fixture for the future PL–LBM sampler.

## Usage

``` r
pl_lbm_toy(n_rankings = 30L, rank_length = 4L, seed = 2027L)
```

## Arguments

- n_rankings:

  Number of ranking rows.

- rank_length:

  Number of reported positions.

- seed:

  Integer seed.

## Value

A `btsbm_ranking_data` object with a `"truth"` attribute containing the
generating ranking and item partitions.
