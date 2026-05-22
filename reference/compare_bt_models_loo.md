# Compare BT models with Pareto-smoothed importance sampling LOO

Convenience wrapper around loo for comparing two log-likelihood matrices
(simple vs clustered BT).

## Usage

``` r
compare_bt_models_loo(simple_llo, cluster_llo)

compare_bt_models_loo(simple_llo, cluster_llo)
```

## Arguments

- simple_llo:

  list returned by
  [`make_bt_simple_loo()`](https://laposanti.github.io/BTSBM/reference/make_bt_simple_loo.md).

- cluster_llo:

  list returned by
  [`make_bt_cluster_loo()`](https://laposanti.github.io/BTSBM/reference/make_bt_cluster_loo.md).

## Value

A list with:

- `simple` — `loo` object for the simple BT.

- `cluster` — `loo` object for the clustered BT–SBM.

- `comparison` — result of `loo::compare_models()`.

A list with:

- `simple` — `loo` object for the simple BT.

- `cluster` — `loo` object for the clustered BT–SBM.

- `comparison` — result of
  [`loo::loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html).
