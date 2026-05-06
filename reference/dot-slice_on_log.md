# Slice sampler on log-scale (internal)

Univariate slice sampling on \\\log a\\. Used internally for
hyperparameter updates.

## Usage

``` r
.slice_on_log(logpost, loga0, w = 1, m = 20L, lower = -Inf, upper = Inf)
```

## Arguments

- logpost:

  Function taking `loga` and returning log-posterior up to a constant.

- loga0:

  Numeric; current log value.

- w:

  Step-out width (default 1).

- m:

  Max step-out steps (default 20).

- lower, upper:

  Hard bounds on `loga`.

## Value

A new `loga` draw.
