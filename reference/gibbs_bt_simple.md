# Simple Bradley–Terry Gibbs sampler (no clustering)

Simple Bradley–Terry Gibbs sampler (no clustering)

## Usage

``` r
gibbs_bt_simple(
  w_ij,
  a = 0.01,
  b = 0.1,
  T_iter = 5000,
  T_burn = 1000,
  verbose = TRUE
)
```

## Arguments

- w_ij:

  integer/numeric \\n \times n\\ wins from i over j (diag = 0,
  nonnegative).

- a, b:

  numeric(1) Gamma(a,b) prior on each \\\lambda_i\\.

- T_iter, T_burn:

  integers; total iterations and burn-in. Require `T_burn < T_iter`.

- verbose:

  logical; print progress every 1000 iterations.

## Value

A list with `lambda_samples` (matrix of size \\S \times n\\,
\\S=T\_{\text{iter}}-T\_{\text{burn}}\\).

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
n <- 6L
w <- matrix(0L, n, n)
w[upper.tri(w)] <- rpois(sum(upper.tri(w)), 2)
w <- w + t(w) - diag(diag(w))
fit <- gibbs_bt_simple(w, a = 1, b = 1, T_iter = 500, T_burn = 100, verbose = FALSE)
} # }
```
