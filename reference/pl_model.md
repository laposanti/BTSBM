# Define a simple Plackett–Luce model

Define a simple Plackett–Luce model

## Usage

``` r
pl_model(
  latent_strength = .default_latent_strength(),
  identifiability = pl_identifiability(),
  ability_prior = NULL
)
```

## Arguments

- latent_strength:

  Gamma prior for positive item strengths.

- identifiability:

  PL reporting constraint.

- ability_prior:

  Deprecated alias for `latent_strength`.

## Value

A `btsbm_model`.
