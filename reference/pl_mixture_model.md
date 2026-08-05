# Define a Plackett–Luce ranking-mixture model

Ranking rows are assigned to latent preference groups. Each group has a
separate vector of item strengths.

## Usage

``` r
pl_mixture_model(
  latent_strength = .default_latent_strength(),
  ranking_clustering = gnedin_prior(gnedin_hyperparameter = 0.5),
  identifiability = pl_identifiability(),
  ability_prior = NULL,
  ranking_prior = NULL
)
```

## Arguments

- latent_strength:

  Gamma prior for positive component strengths.

- ranking_clustering:

  Prior on latent ranking groups.

- identifiability:

  PL reporting constraint.

- ability_prior:

  Deprecated alias for `latent_strength`.

- ranking_prior:

  Deprecated alias for `ranking_clustering`.

## Value

A `btsbm_model`.
