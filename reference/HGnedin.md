# Gnedin prior normalization weight K(n,k)

Computes the unnormalized mass function term used in Gnedin-type priors.

## Usage

``` r
HGnedin(n, k, gamma = 0.5)
```

## Arguments

- n:

  integer(1) total items.

- k:

  integer vector of block counts.

- gamma:

  numeric(1) parameter \\\gamma \> 0\\.

## Value

Numeric vector of weights \\K(n,k)\\.

## References

Legramanti, S., Rigon, T., Durante, D., Dunson, D.B., 2022. Extended
stochastic block models with application to criminal networks. The
Annals of Applied Statistics 16. https://doi.org/10.1214/21-aoas1595
