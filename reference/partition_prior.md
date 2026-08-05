# Legacy partition-prior interface

Legacy partition-prior interface

## Usage

``` r
partition_prior(
  type = c("gnedin", "dp", "py", "dm"),
  gamma = 0.5,
  alpha = 1,
  sigma = 0.2,
  beta = 1,
  K = NULL
)
```

## Arguments

- type:

  One of `"gnedin"`, `"dp"`, `"py"`, or `"dm"`.

- gamma:

  Gnedin hyperparameter.

- alpha:

  DP/Pitman–Yor concentration.

- sigma:

  Pitman–Yor discount.

- beta:

  Finite-Dirichlet concentration.

- K:

  Maximum finite-Dirichlet groups.

## Value

An object of class `btsbm_partition_prior`.
