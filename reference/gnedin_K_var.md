# Variance of the number of clusters under the Gnedin prior

Compute the closed-form variance of \\K_n\\ under the Gnedin finite-type
model (\\\sigma=-1\\).

## Usage

``` r
gnedin_K_var(n, gamma)
```

## Arguments

- n:

  Integer vector of sample sizes (each \\\ge 1\\).

- gamma:

  Numeric vector of Gnedin parameters, each in \\(0,1)\\.

## Value

A numeric vector with \\\mathrm{Var}(K_n)\\ for each pair `(n, gamma)`.

## Details

Using the factorial–ordinary moment relation
\\\mathbb{E}\[K_n^2\]=\mathbb{E}\[K_n(K_n-1)\]+\mathbb{E}\[K_n\]\\, the
variance is \$\$ \mathrm{Var}(K_n) \\=\\ \mathbb{E}\[K_n(K_n-1)\] +
\mathbb{E}\[K_n\] - \\\mathbb{E}\[K_n\]\\^2, \$\$ where \$\$
\mathbb{E}\[K_n(K_n-1)\] \\=\\
n(n-1)(1-\gamma)\\\frac{\Gamma(n)\\\Gamma(1+\gamma)}{\Gamma(n+\gamma)}.
\$\$ The implementation uses `lgamma` for numerical stability and is
vectorized over `n` and `gamma` (with standard R recycling rules).

## References

Gnedin, A. (2010). A species sampling model with finitely many types.
*Electronic Communications in Probability*, 15, 79–88.

Favaro, S., Lijoi, A., & Prünster, I. (2013). Extending the class of
Gibbs-type priors: theoretical properties and new examples. *Annals of
Applied Probability*, 23(4), 1729–1754.

Pitman, J. (2006). *Combinatorial Stochastic Processes*. Springer.

## See also

[`gnedin_K_mean()`](https://laposanti.github.io/BTSBM/reference/gnedin_K_mean.md)

## Examples

``` r
# Scalar inputs
gnedin_K_var(105, 0.8)
#> [1] 45.95132

# Consistency check: variance is nonnegative
all(gnedin_K_var(20, c(0.3, 0.5, 0.8)) >= 0)
#> [1] TRUE
```
