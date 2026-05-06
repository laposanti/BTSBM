# Format a block-count posterior distribution as a pretty- table

Given the output of relabel_by_lambda, this function produces a pretty
table with the posterior distribution of K, the number of blocks induced
by the partition x. It formats the result with
[`knitr::kable`](https://rdrr.io/pkg/knitr/man/kable.html) and
[`kableExtra::kable_styling`](https://rdrr.io/pkg/kableExtra/man/kable_styling.html).

## Usage

``` r
pretty_table_K_distribution(
  post,
  format = "html",
  p_min = 0.05,
  digits = 3,
  latex_options = c("striped", "hold_position")
)
```

## Arguments

- post:

  The output of relabel_by_lambda()

- format:

  Character string. The format of the output table. Default is "html".
  Switch to "latex" for a document-ready table.

- p_min:

  Numeric scalar. Probabilities \\\le\\ this threshold are dropped
  before pivoting. Defaults to `0.05`.

- digits:

  Integer passed to
  [`knitr::kable`](https://rdrr.io/pkg/knitr/man/kable.html) for
  rounding. Defaults to `3`.

- latex_options:

  Character vector of LaTeX options passed to
  [`kableExtra::kable_styling`](https://rdrr.io/pkg/kableExtra/man/kable_styling.html).
  Defaults to `c("striped","hold_position")`.

## Value

A `kableExtra` / `knitr_kable` object (suitable for LaTeX or HTML,
depending on `format`).

## Details

The function produces a pretty table with the posterior distribution of
K, the number of blocks induced by the partition x. it filters the
clusters that have a p(K \> `p_min`).
