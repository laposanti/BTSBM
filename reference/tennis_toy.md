# Tennis-league toy paired-comparison data

A synthetic, tiered round-robin tournament for demonstrating the
Bradley–Terry and BT–SBM interfaces. It contains no real player results.

## Usage

``` r
tennis_toy(matches_per_pair = 8L, seed = 2031L)
```

## Arguments

- matches_per_pair:

  Positive number of matches played by each pair.

- seed:

  Integer seed.

## Value

A btsbm_pairwise_data object. Its truth attribute contains the
generating player strengths and latent tiers.
