# Define a simple Bradley–Terry model

Define a simple Bradley–Terry model

## Usage

``` r
bt_model(
  latent_strength = .default_latent_strength(shape = 0.01, rate = 0.1),
  ability_prior = NULL
)
```

## Arguments

- latent_strength:

  Gamma prior for positive item strengths.

- ability_prior:

  Deprecated alias for `latent_strength`.

## Value

A `btsbm_model`.
