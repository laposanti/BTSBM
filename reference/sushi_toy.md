# Sushi-inspired toy ranking data

A reproducible, synthetic teaching dataset with two preference profiles
over ten sushi types. It is designed to demonstrate
[`as_rankings()`](https://laposanti.github.io/BTSBM/reference/as_rankings.md),
[`pl_model()`](https://laposanti.github.io/BTSBM/reference/pl_model.md),
PL mixtures, and PL–LBM data shapes without redistributing
respondent-level data from the Sushi Preference Dataset.

## Usage

``` r
sushi_toy(n_rankings = 40L, rank_length = 5L, seed = 2026L)
```

## Arguments

- n_rankings:

  Number of synthetic rankings.

- rank_length:

  Number of reported positions per ranking.

- seed:

  Integer seed.

## Value

A `btsbm_ranking_data` object. Its `"truth"` attribute contains
synthetic ranking-cluster labels and the generating ability matrix.

## References

Kamishima, T. (2003). Nantonac collaborative filtering: Recommendation
based on order responses. *Proceedings of the Ninth ACM SIGKDD
International Conference on Knowledge Discovery and Data Mining*. The
original Sushi Preference Dataset must be obtained from its authors; its
licence does not allow redistribution.
