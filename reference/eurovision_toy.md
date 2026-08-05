# Eurovision-inspired toy ranking data

A synthetic preference-ranking dataset for a song-contest setting. Two
latent voting profiles produce different rankings of eight fictional
acts. It is designed to demonstrate
[`pl_mixture_model()`](https://laposanti.github.io/BTSBM/reference/pl_mixture_model.md)
without redistributing data from a real contest.

## Usage

``` r
eurovision_toy(n_rankings = 40L, rank_length = 5L, seed = 2030L)
```

## Arguments

- n_rankings:

  Number of synthetic jury or viewer rankings.

- rank_length:

  Number of reported positions per ranking.

- seed:

  Integer seed.

## Value

A btsbm_ranking_data object. Its truth attribute contains the synthetic
voter-group labels and generating latent strengths.
