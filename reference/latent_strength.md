# Define a Gamma prior for positive latent strengths

A Plackett–Luce likelihood only identifies relative strengths. This
prior is defined on the positive, unnormalised strengths used by the
sampler;
[`pl_identifiability()`](https://laposanti.github.io/BTSBM/reference/pl_identifiability.md)
controls the scale used when strengths are reported.

## Usage

``` r
latent_strength(shape = 1, rate = 1)

gamma_ability(shape = 1, rate = 1)
```

## Arguments

- shape:

  Positive Gamma shape for each latent strength.

- rate:

  Positive Gamma rate for each latent strength.

## Value

An object of class `btsbm_strength_prior`.
