# Expected number of clusters under the Gnedin prior

Compute the closed-form mean of \\K_n\\ (number of clusters) under the
Gnedin finite-type model (\\\sigma=-1\\).

## Usage

``` r
gnedin_K_mean(n, gamma)
```

## Arguments

- n:

  Integer vector of sample sizes (each \\\ge 1\\).

- gamma:

  Numeric vector of Gnedin parameters, each in \\(0,1)\\.

## Value

A numeric vector with \\\mathbb{E}\[K_n\]\\ for each pair `(n, gamma)`.

## Details

For sample size \\n \ge 1\\ and parameter \\\gamma \in (0,1)\\, the mean
is \$\$ \mathbb{E}\[K_n\] \\=\\
\frac{\Gamma(n+1)\\\Gamma(1+\gamma)}{\Gamma(n+\gamma)}. \$\$ The
implementation uses `lgamma` for numerical stability and is vectorized
over `n` and `gamma` (with standard R recycling rules).

The formula follows directly from standard Gibbs–type manipulations
using factorial moments and the Chu–Vandermonde identity specialized to
\\\sigma=-1\\.

## References

Gnedin, A. (2010). A species sampling model with finitely many types.
*Electronic Communications in Probability*, 15, 79–88.

Favaro, S., Lijoi, A., & Prünster, I. (2013). Extending the class of
Gibbs-type priors: theoretical properties and new examples. *Annals of
Applied Probability*, 23(4), 1729–1754.

Pitman, J. (2006). *Combinatorial Stochastic Processes*. Springer.

## See also

[`gnedin_K_var()`](https://laposanti.github.io/BTSBM/reference/gnedin_K_var.md)

## Examples

``` r
# Scalar inputs
gnedin_K_mean(105, 0.8)
#> [1] 2.36427

# Vectorized over gamma
gnedin_K_mean(50, c(0.3, 0.5, 0.8))
#> [1] 13.906329  6.282256  2.039934

# Vectorized over n and gamma (recycling rules)
gnedin_K_mean(c(20, 50, 100), 0.5)
#> [1] 3.988173 6.282256 8.873354
```
