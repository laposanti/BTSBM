# Define a Plackett–Luce latent block model

The PL–LBM jointly partitions ranking rows and items. Each pair of
latent groups has its own positive strength.

## Usage

``` r
pl_lbm_model(
  latent_strength = .default_latent_strength(),
  ranking_clustering = gnedin_prior(gnedin_hyperparameter = 0.5),
  item_clustering = gnedin_prior(gnedin_hyperparameter = 0.5),
  identifiability = pl_identifiability(),
  ability_prior = NULL,
  ranking_prior = NULL,
  item_prior = NULL
)
```

## Arguments

- latent_strength:

  Gamma prior for positive cell strengths.

- ranking_clustering:

  Prior on latent ranking groups.

- item_clustering:

  Prior on latent item blocks.

- identifiability:

  PL reporting constraint.

- ability_prior:

  Deprecated alias for `latent_strength`.

- ranking_prior:

  Deprecated alias for `ranking_clustering`.

- item_prior:

  Deprecated alias for `item_clustering`.

## Value

A `btsbm_model`.
