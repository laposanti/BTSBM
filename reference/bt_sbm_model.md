# Define a Bradley–Terry stochastic block model

Define a Bradley–Terry stochastic block model

## Usage

``` r
bt_sbm_model(
  latent_strength = .default_latent_strength(shape = 4, rate = exp(digamma(4))),
  item_clustering = gnedin_prior(gnedin_hyperparameter = 0.5),
  ability_prior = NULL,
  item_prior = NULL
)
```

## Arguments

- latent_strength:

  Gamma prior for positive block strengths.

- item_clustering:

  Prior on latent item blocks.

- ability_prior:

  Deprecated alias for `latent_strength`.

- item_prior:

  Deprecated alias for `item_clustering`.

## Value

A `btsbm_model`.
