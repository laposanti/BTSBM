# Simulate strict Plackett–Luce rankings from a weight matrix

Simulate strict Plackett–Luce rankings from a weight matrix

## Usage

``` r
.simulate_pl_rankings(strength, rank_length = ncol(strength))
```

## Arguments

- strength:

  Matrix with one positive latent-strength vector per ranking row.

- rank_length:

  Number of retained ranking positions.

## Value

Integer ranking matrix.
