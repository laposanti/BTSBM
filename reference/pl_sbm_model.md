# Define a Plackett–Luce stochastic block model

Define a Plackett–Luce stochastic block model

## Usage

``` r
pl_sbm_model(
  latent_strength = .default_latent_strength(),
  item_clustering = gnedin_prior(gnedin_hyperparameter = 0.5),
  identifiability = pl_identifiability(),
  ability_prior = NULL,
  item_prior = NULL
)
```

## Arguments

- latent_strength:

  Gamma prior for positive block strengths.

- item_clustering:

  Prior on latent item blocks. Use
  [`gnedin_prior()`](https://laposanti.github.io/BTSBM/reference/gnedin_prior.md),
  [`dirichlet_process_prior()`](https://laposanti.github.io/BTSBM/reference/dirichlet_process_prior.md),
  [`pitman_yor_prior()`](https://laposanti.github.io/BTSBM/reference/pitman_yor_prior.md),
  or
  [`finite_partition()`](https://laposanti.github.io/BTSBM/reference/finite_partition.md).

- identifiability:

  PL reporting constraint.

- ability_prior:

  Deprecated alias for `latent_strength`.

- item_prior:

  Deprecated alias for `item_clustering`.

## Value

A `btsbm_model`.
