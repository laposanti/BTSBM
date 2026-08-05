# Choose a prior for a latent partition

This is the common prior selector for item blocks in PL–SBM and ranking
groups in a PL mixture. `family = "gnedin"` is a flexible finite-cluster
prior and is the default. A Dirichlet process and Pitman–Yor process
allow a potentially unbounded number of groups. A finite Dirichlet prior
fixes an upper bound through `max_clusters`.

## Usage

``` r
clustering_prior(
  family = c("gnedin", "dirichlet_process", "pitman_yor", "finite_dirichlet"),
  gnedin_hyperparameter = 0.5,
  concentration = 1,
  discount = 0.2,
  max_clusters = NULL
)
```

## Arguments

- family:

  One of `"gnedin"`, `"dirichlet_process"`, `"pitman_yor"`, or
  `"finite_dirichlet"`.

- gnedin_hyperparameter:

  Gnedin hyperparameter in `(0, 1)`. Larger values favour more occupied
  groups a priori.

- concentration:

  Concentration for the Dirichlet-process, Pitman–Yor, or
  finite-Dirichlet prior.

- discount:

  Pitman–Yor discount in `[0, 1)`.

- max_clusters:

  Required maximum number of groups for `family = "finite_dirichlet"`.

## Value

An object of class `btsbm_partition_prior`.
